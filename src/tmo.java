package tmo;

import com.stata.sfi.Mata;
import com.stata.sfi.Data;
import com.stata.sfi.Matrix;
import com.stata.sfi.SFIToolkit;
import com.stata.sfi.Scalar;
import java.util.BitSet;
import java.util.Arrays;
import java.lang.management.ManagementFactory;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Thresholded anchor variance paths and target-pair diagnostics.
 *
 * For every outcome a in 0..D, treat a as the anchor: orthogonalize the other
 * outcomes on it, screen location pairs on the resulting statistic, and
 * accumulate that anchor's thresholded variance over a grid of thresholds.
 *
 * Pair screens use a rank-two identity. With Xs the n-by-p screening
 * residuals and h_a = gram[pool,a]/ss[a],
 *
 *     sum_j eps~(j;a)_i eps~(j;a)_j  =  [XX' - z_a x_a' - x_a z_a']_{i,j}
 *     where z_a = Xs[,pool] h_a - x_a (h_a'h_a)/2 and XX' uses pool only.
 * The default pool contains all outcomes; an anchor's own projected residual
 * is zero. Screening inputs carry the requested column scaling, before
 * anchor projection.
 *
 * All p anchors share one Gram matrix and one strict triangle of pair dot
 * products. The latter is reused between endpoint calibration and paths.
 * Repeated calls reuse work arrays and daemon workers. Explicit simulation
 * contexts may update the triangle across overlapping outcome windows;
 * uncertain rank/threshold comparisons always recompute the original dot.
 *
 * Inputs are read as GLOBAL Mata matrices in one bulk call each:
 *   __tmo_Xs     n x p  screening residuals
 *   __tmo_W      n x p  weights, already multiplied by xtilde
 *   __tmo_dgrid  G x 1  ascending positive thresholds (explicit-threshold mode)
 * plus Stata scalars __tmo_corrscale (1 = correlation screening) and
 * __tmo_threads (0 or missing = all available cores). Optional flags, default 0:
 *   __tmo_auxonly=1  restrict the screening pool to auxiliary columns 1..d
 *   __tmo_masked=1   protect the strict-lower pairs in __tmo_Bpairs (q x 2,
 *                    one-based indices) and calibrate on eligible off-diagonal
 *                    pairs only. Unmasked endpoint calibration uses all n^2 entries.
 *   __tmo_Bgroup     optional n x 1 contiguous positive group codes: a compact
 *                    single partition, requiring an empty __tmo_Bpairs.
 *   __tmo_need_counts=1 returns target screened counts (1 x G) and the number
 *                    of unsupported unordered target pairs, including protected.
 * Kernel API is 4, mask API is 2; counts API 1 marks requested count outputs.
 * Diagnostic API 2 prepares one immutable snapshot, streams pair blocks, and
 * releases it explicitly; a new fit also invalidates any diagnostic snapshot.
 *   __tmo_endpoint=1 build delta_k=sqrt(2*eta_k*v0hat*logfactor) from the
 *                    G x 1 grid __tmo_etagrid and __tmo_logfactor>0. This mode
 *                    requires correlation screens and computes anchor-0 IQR
 *                    from eligible pairs if masked, otherwise all n^2 entries.
 * Endpoint outputs: scalars __tmo_v0hat, __tmo_q25, __tmo_q75; the G x 1
 * Stata matrix __tmo_dgrid_out. Its quantiles use the inverse empirical CDF.
 *
 * Output is the Stata matrix __tmo_paths, p x G, each entry
 *     sum_i W[i,a]^2 + protected pair products + eligible thresholded products
 * i.e. Vhat's numerator before dividing by (xtilde'xtilde)^2, with the
 * diagonal always retained.
 *
 * The reduction is deterministic: each thread accumulates its own bins and the
 * per-thread arrays are summed back in fixed index order. Runs with the same
 * thread count are deterministic; changing thread counts may change last bits.
 */
public final class tmo {

    private tmo() { }

    private static final Workspace workspace = new Workspace();
    private static ExecutorService workerPool;
    private static int workerPoolSize;
    private static boolean reuseWorkers = true;
    private static DiagnosticState diagnosticState;
    private static long activeRunStart;

    private static final class DiagnosticState {
        int n, p, first;
        boolean explicit;
        double[] x, residual, z, rsd, gram, ss;
        boolean[] observed;
    }

    /** One bounded workspace per Stata JVM; entrypoints are serialized. */
    private static final class Workspace {
        double[] gram, ss, hh, z, gd, rsd, lower, common, wb, wbError, missingResidual;
        double[][] partial, partialError, missingWorkers;
        boolean[] observed;
        double[] previousX;
        int previousN, previousP, shifts;
        double context, window, maximumX;
        boolean rollingValid;
        long hits, rebuilds;

        void invalidate() { rollingValid = false; previousX = null; shifts = 0; }
        void release() {
            gram = ss = hh = z = gd = rsd = lower = common = wb = wbError = missingResidual = null;
            partial = partialError = missingWorkers = null;
            observed = null;
            invalidate();
        }
        long retainedBytes() {
            long count = 0;
            for (double[] a : new double[][] {gram, ss, hh, z, gd, rsd,
                    lower, common, wb, wbError, previousX, missingResidual}) if (a != null) count += a.length;
            if (partial != null) for (double[] a : partial) if (a != null) count += a.length;
            if (partialError != null) for (double[] a : partialError) if (a != null) count += a.length;
            if (missingWorkers != null) for (double[] a : missingWorkers) if (a != null) count += a.length;
            return 8L * count + (observed == null ? 0 : observed.length);
        }
    }

    private static double[] array(double[] old, int size, boolean clear) {
        if (old == null || old.length != size) return new double[size];
        if (clear) Arrays.fill(old, 0.0);
        return old;
    }

    private static long reusable(double[] value, long length) {
        return value != null && value.length == length ? 8L * length : 0L;
    }

    private static long reusableBytes(int n, int p, int G, int threads, boolean endpoint,
            boolean missing, boolean cachePairs, long eligible, int missingWorkers) {
        if (!reuseWorkers) return 0;
        final long cells = (long) n * p, bins = (long) p * (G + 1);
        long bytes = endpoint ? reusable(workspace.lower, eligible) : 0;
        bytes += reusable(workspace.wb, bins) + reusable(workspace.wbError, bins);
        if (missing) {
            if (workspace.observed != null && workspace.observed.length == cells) bytes += cells;
            bytes += reusable(workspace.missingResidual, cells);
            if (workspace.missingWorkers != null && workspace.missingWorkers.length == missingWorkers)
                for (double[] buffer : workspace.missingWorkers) bytes += reusable(buffer, cells);
        } else {
            bytes += reusable(workspace.gram, (long) p * p) + reusable(workspace.ss, p)
                + reusable(workspace.hh, p) + reusable(workspace.z, cells)
                + reusable(workspace.gd, n) + reusable(workspace.rsd, cells);
            if (workspace.partial != null && workspace.partial.length == threads)
                for (double[] buffer : workspace.partial) bytes += reusable(buffer, bins);
            if (workspace.partialError != null && workspace.partialError.length == threads)
                for (double[] buffer : workspace.partialError) bytes += reusable(buffer, bins);
            if (cachePairs) bytes += reusable(workspace.common, (long) n * (n - 1) / 2)
                + reusable(workspace.previousX, cells);
        }
        return bytes;
    }

    public static synchronized int run(String[] args) {
        final long start = System.nanoTime();
        activeRunStart = start;
        try {
            diagnosticState = null;
            Scalar.setValue("__tmo_diagnostic_ready", 0.0);
            Scalar.setValue("__tmo_diagnostic_preparations", 0.0);
            Scalar.setValue("__tmo_kernel_api", 4.0);
            Scalar.setValue("__tmo_explicit_screen", 0.0);
            Scalar.setValue("__tmo_java_memory_failure", 0.0);
            Scalar.setValue("__tmo_mask_api", 2.0);
            Scalar.setValue("__tmo_counts_api", 0.0);
            Scalar.setValue("__tmo_calibration_error", 0.0);
            Scalar.setValue("__tmo_preparation_ms", 0.0);
            Scalar.setValue("__tmo_calibration_ms", 0.0);
            Scalar.setValue("__tmo_paths_ms", 0.0);
            return compute();
        } catch (OutOfMemoryError | MemoryLimitException t) {
            workspace.release();
            Scalar.setValue("__tmo_java_memory_failure", 1.0);
            SFIToolkit.errorln("tmo: insufficient Java memory for exact screening. "
                + t.getMessage() + " Reduce the location count or increase the Java heap; "
                + "the larger Mata fallback must not be attempted.");
            return error(909);
        } catch (ProcedureException t) {
            workspace.invalidate();
            SFIToolkit.errorln("tmo: " + t.getMessage());
            return error(459);
        } catch (Throwable t) {
            workspace.invalidate();
            SFIToolkit.errorln("tmo: " + t);
            return error(498);
        } finally {
            Scalar.setValue("__tmo_compute_ms", (System.nanoTime() - start) / 1e6);
            Scalar.setValue("__tmo_workspace_bytes", workspace.retainedBytes());
            Scalar.setValue("__tmo_heap_used_bytes", Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
            Scalar.setValue("__tmo_pool_workers", workerPoolSize);
        }
    }

    private static int error(int rc) {
        return rc;
    }

    /** Best-effort host-RAM preflight; this is not a memory reservation. */
    public static synchronized int hostMemory(String[] args) {
        Scalar.setValue("__tmo_host_free_bytes", Data.getMissingValue());
        try {
            final Object bean = ManagementFactory.getOperatingSystemMXBean();
            final Class<?> extended = Class.forName("com.sun.management.OperatingSystemMXBean");
            if (extended.isInstance(bean)) {
                final long bytes = ((Number) extended.getMethod("getFreePhysicalMemorySize")
                    .invoke(bean)).longValue();
                if (bytes > 0) Scalar.setValue("__tmo_host_free_bytes", (double) bytes);
            }
        } catch (Throwable unavailable) { /* Optional on non-HotSpot platforms. */ }
        return 0;
    }

    /** Target diagnostics for one-based pairs supplied as Mata column vectors. */
    public static synchronized int targetScores(String[] args) {
        try {
            Scalar.setValue("__tmo_java_memory_failure", 0.0);
            return computeTargetScores();
        } catch (OutOfMemoryError | MemoryLimitException failure) {
            Scalar.setValue("__tmo_java_memory_failure", 1.0);
            SFIToolkit.errorln("tmo: insufficient Java memory for pair diagnostics. "
                + failure.getMessage());
            return error(909);
        } catch (Throwable failure) {
            SFIToolkit.errorln("tmo targetScores: " + failure);
            return error(498);
        }
    }

    /** A single export owns this immutable snapshot until release or the next fit. */
    public static synchronized int prepareTargetScores(String[] args) {
        diagnosticState = null;
        Scalar.setValue("__tmo_diagnostic_ready", 0.0);
        Scalar.setValue("__tmo_diagnostic_api", 2.0);
        try {
            diagnosticState = prepareDiagnosticState();
            Scalar.setValue("__tmo_diagnostic_ready", 1.0);
            return 0;
        } catch (OutOfMemoryError | MemoryLimitException failure) {
            Scalar.setValue("__tmo_java_memory_failure", 1.0);
            SFIToolkit.errorln("tmo: insufficient Java memory for target diagnostics");
            return 909;
        } catch (Throwable failure) {
            SFIToolkit.errorln("tmo target preparation: " + failure);
            return 498;
        }
    }

    public static synchronized int releaseTargetScores(String[] args) {
        diagnosticState = null;
        Scalar.setValue("__tmo_diagnostic_ready", 0.0);
        return 0;
    }

    private static DiagnosticState prepareDiagnosticState() {
        final DiagnosticState state = new DiagnosticState();
        state.n = (int) Mata.getMataRowTotal("__tmo_Xs");
        state.p = (int) Mata.getMataColTotal("__tmo_Xs");
        final int n = state.n, p = state.p;
        if (n < 2 || p < 2) throw new IllegalArgumentException("bad screening dimensions");
        state.first = optionalScalar("__tmo_auxonly", 0.0) == 1.0 ? 1 : 0;
        state.explicit = optionalScalar("__tmo_hasmiss", 0.0) == 1.0
            || optionalScalar("__tmo_explicit_screen", 0.0) == 1.0;
        checkMemory(n, p, 1, 1, false, state.explicit);
        state.x = Mata.getMataReal("__tmo_Xs");
        final double[] x = state.x;
        if (!state.explicit) {
            state.gram = new double[p]; state.ss = new double[1];
            state.z = new double[n]; state.rsd = new double[n];
            for (int i = 0; i < n; i++) {
                final int base = i * p;
                for (int b = 0; b < p; b++) state.gram[b] += x[base] * x[base + b];
            }
            state.ss[0] = state.gram[0];
            double hh = 0.0;
            for (int b = state.first; b < p; b++) {
                final double h = state.gram[b] / state.ss[0]; hh += h * h;
            }
            for (int i = 0; i < n; i++) {
                final int base = i * p;
                double acc = 0.0, gd = 0.0;
                for (int b = state.first; b < p; b++) {
                    acc += x[base + b] * (state.gram[b] * (1.0 / state.ss[0]));
                    gd += x[base + b] * x[base + b];
                }
                state.z[i] = acc - x[base] * hh * .5;
                final double subtract = 2.0 * state.z[i] * x[base], v = (gd - subtract) / (p - 1);
                if (illConditionedNorm(gd, subtract, p - 1)) state.explicit = true;
                state.rsd[i] = v > 0.0 ? 1.0 / Math.sqrt(v) : 0.0;
            }
        }
        if (state.explicit) {
            state.observed = new boolean[x.length]; state.residual = new double[x.length];
            for (int k = 0; k < x.length; k++) {
                state.observed[k] = validNumber(x[k]);
                if (!state.observed[k]) x[k] = 0.0;
            }
            missingResiduals(x, state.observed, state.residual, n, p, 0, state.first);
        }
        Scalar.setValue("__tmo_diagnostic_api", 2.0);
        Scalar.setValue("__tmo_diagnostic_preparations", optionalScalar("__tmo_diagnostic_preparations", 0.0) + 1.0);
        return state;
    }

    private static int computeTargetScores() {
        final DiagnosticState state = diagnosticState != null ? diagnosticState : prepareDiagnosticState();
        final int n = state.n, p = state.p, first = state.first;
        final boolean missing = state.explicit;
        final double[] x = state.x, residual = state.residual, z = state.z, rsd = state.rsd;
        final boolean[] observed = state.observed;
        if (Mata.getMataRowTotal("__tmo_Xs") != n || Mata.getMataColTotal("__tmo_Xs") != p)
            throw new IllegalArgumentException("target diagnostic snapshot dimensions changed; prepare it again");
        final double[] pairI = Mata.getMataReal("__tmo_pair_i");
        final double[] pairJ = Mata.getMataReal("__tmo_pair_j");
        if (pairI.length != pairJ.length) throw new IllegalArgumentException("pair vectors differ in length");
        final int count = pairI.length;
        execute("mata: __tmo_pair_scores=J(" + count + ",1,.); "
            + "__tmo_pair_supported=J(" + count + ",1,.); "
            + "__tmo_pair_shared=J(" + count + ",1,.)");
        final int[] details = new int[2];
        try {
            for (int lo = 0; lo < count; lo += 256) {
                final int rows = Math.min(256, count - lo);
                if (Matrix.createMatrix("__tmo_pair_chunk", rows, 3, 0.0) != 0)
                    throw new IllegalStateException("cannot allocate the diagnostic bridge");
                for (int k = 0; k < rows; k++) {
                    final double ri = pairI[lo + k], rj = pairJ[lo + k];
                    if (ri != Math.rint(ri) || rj != Math.rint(rj)
                            || ri < 1 || ri > n || rj < 1 || rj > n)
                        throw new IllegalArgumentException("pair index outside the screening sample");
                    // Match the variance loop's strict lower-triangle arithmetic.
                    final int i = (int) Math.max(ri, rj) - 1, j = (int) Math.min(ri, rj) - 1;
                    final int bi = i * p, bj = j * p;
                    double score = 1.0;
                    int shared = p - 1;
                    boolean supported = true;
                    if (missing) {
                        final double candidate = missingScore(residual, observed, i, j, p, 0, first, details);
                        shared = details[0];
                        if (i != j) {
                            supported = details[1] == 1;
                            score = candidate;
                        }
                    } else if (i != j) {
                        double common = 0.0;
                        for (int b = first; b < p; b++) common += x[bi + b] * x[bj + b];
                        score = (common - z[i] * x[bj] - x[bi] * z[j]) * (1.0 / (p - 1));
                        if (illConditionedPair(common, z[i] * x[bj], x[bi] * z[j])
                                || nearPerfect(score * rsd[i] * rsd[j]))
                            score = explicitCompletePair(x, state.gram, state.ss, bi, bj, p, 0, first);
                        else score = boundCorrelation(score * rsd[i] * rsd[j]);
                    }
                    Matrix.storeMatrixAt("__tmo_pair_chunk", k, 0, score);
                    Matrix.storeMatrixAt("__tmo_pair_chunk", k, 1, supported ? 1.0 : 0.0);
                    Matrix.storeMatrixAt("__tmo_pair_chunk", k, 2, shared);
                }
                final String range = "[|" + (lo + 1) + ",1 \\ " + (lo + rows) + ",1|]";
                execute("mata: __tmo_pair_scores" + range + "=st_matrix(\"__tmo_pair_chunk\")[,1]; "
                    + "__tmo_pair_supported" + range + "=st_matrix(\"__tmo_pair_chunk\")[,2]; "
                    + "__tmo_pair_shared" + range + "=st_matrix(\"__tmo_pair_chunk\")[,3]");
            }
        } finally {
            SFIToolkit.executeCommand("capture matrix drop __tmo_pair_chunk", false);
        }
        return 0;
    }

    private static void execute(String command) {
        final int rc = SFIToolkit.executeCommand(command, false);
        if (rc != 0) throw new IllegalStateException("diagnostic bridge returned " + rc);
    }

    private static void createMatrix(String name, int rows, int columns) {
        final int rc = Matrix.createMatrix(name, rows, columns, 0.0);
        if (rc != 0) throw new ProcedureException("cannot allocate output matrix " + name + " (return code " + rc + ")");
    }

    private static final class MemoryLimitException extends RuntimeException {
        MemoryLimitException(String message) { super(message); }
    }

    private static final class ProcedureException extends RuntimeException {
        ProcedureException(String message) { super(message); }
    }

    private static int pairIndex(int i, int j) {
        return (int) ((long) i * (i - 1) / 2 + j);
    }

    private static final class ProtectedPairs {
        final boolean enabled;
        final BitSet bits;
        final long[] eligibleOffsets;
        final int[] groups, members, groupOffsets;
        final long eligibleCount;

        ProtectedPairs(int n, boolean enabled) {
            this.enabled = enabled;
            final long pairs = (long) n * (n - 1) / 2;
            if (!enabled) {
                bits = null;
                groups = members = groupOffsets = null;
                eligibleOffsets = null;
                eligibleCount = pairs;
                return;
            }
            eligibleOffsets = new long[n + 1];
            final long groupRows = Mata.getMataRowTotal("__tmo_Bgroup");
            if (groupRows > 0) {
                if (groupRows != n || Mata.getMataColTotal("__tmo_Bgroup") != 1
                        || Mata.getMataRowTotal("__tmo_Bpairs") != 0)
                    throw new ProcedureException("__tmo_Bgroup must be n by 1 with an empty __tmo_Bpairs");
                bits = null;
                groups = new int[n];
                final int[] sizes = new int[n];
                final double[] suppliedGroups = Mata.getMataReal("__tmo_Bgroup");
                int count = 0;
                for (int i = 0; i < n; i++) {
                    final double code = suppliedGroups[i];
                    if (!validNumber(code) || code != Math.rint(code) || code < 1 || code > n)
                        throw new ProcedureException("protected group codes must be contiguous positive integers");
                    final int group = (int) code - 1;
                    groups[i] = group;
                    count = Math.max(count, group + 1);
                    eligibleOffsets[i + 1] = eligibleOffsets[i] + i - sizes[group]++;
                }
                groupOffsets = new int[count + 1];
                for (int group = 0; group < count; group++) {
                    if (sizes[group] == 0)
                        throw new ProcedureException("protected group codes must be contiguous positive integers");
                    groupOffsets[group + 1] = groupOffsets[group] + sizes[group];
                }
                members = new int[n];
                Arrays.fill(sizes, 0);
                for (int i = 0; i < n; i++) {
                    final int group = groups[i];
                    members[groupOffsets[group] + sizes[group]++] = i;
                }
                eligibleCount = eligibleOffsets[n];
                return;
            }
            groups = members = groupOffsets = null;
            bits = new BitSet();
            final double[] supplied = Mata.getMataRowTotal("__tmo_Bpairs") == 0
                ? new double[0] : Mata.getMataReal("__tmo_Bpairs");
            final int[] excluded = new int[n];
            for (int k = 0; k < supplied.length; k += 2) {
                final double ri = supplied[k], rj = supplied[k + 1];
                if (!validNumber(ri) || !validNumber(rj) || ri != Math.rint(ri)
                        || rj != Math.rint(rj) || ri <= rj || rj < 1 || ri > n)
                    throw new ProcedureException("protected pairs must be one-based indices with 1 <= j < i <= n");
                final int i = (int) ri - 1, j = (int) rj - 1;
                final int index = pairIndex(i, j);
                if (bits.get(index)) throw new ProcedureException("protected pairs contain a duplicate");
                bits.set(index);
                excluded[i]++;
            }
            for (int i = 0; i < n; i++) eligibleOffsets[i + 1] = eligibleOffsets[i] + i - excluded[i];
            eligibleCount = eligibleOffsets[n];
        }

        boolean contains(int i, int j) {
            return enabled && (groups != null ? groups[i] == groups[j] : bits.get(pairIndex(i, j)));
        }
        int offset(int i) { return (int) (enabled ? eligibleOffsets[i] : (long) i * (i - 1) / 2); }
    }

    /** Neumaier sums plus FMA product residuals preserve cancellation terms. */
    private static void add(double[] sum, double[] correction, int k, double value) {
        final double old = sum[k], next = old + value;
        if (!Double.isFinite(next)) throw new ProcedureException("variance accumulation exceeds the finite numeric range");
        correction[k] += Math.abs(old) >= Math.abs(value) ? (old - next) + value : (value - next) + old;
        sum[k] = next;
    }

    private static void addProduct(double[] sum, double[] correction, int k,
            double left, double right, boolean twice) {
        final double product = left * right;
        if (!Double.isFinite(product)) throw new ProcedureException("variance product exceeds the finite numeric range");
        final double residual = Math.fma(left, right, -product);
        add(sum, correction, k, product);
        add(sum, correction, k, residual);
        if (twice) {
            add(sum, correction, k, product);
            add(sum, correction, k, residual);
        }
    }

    private static double[][] baseline(double[] w, int n, int p, ProtectedPairs mask) {
        final double[] base = new double[p], error = new double[p];
        if (mask.groups != null) {
            final double[] scores = new double[p], scoreError = new double[p];
            for (int group = 0; group + 1 < mask.groupOffsets.length; group++) {
                Arrays.fill(scores, 0.0);
                Arrays.fill(scoreError, 0.0);
                for (int k = mask.groupOffsets[group]; k < mask.groupOffsets[group + 1]; k++) {
                    final int row = mask.members[k] * p;
                    for (int a = 0; a < p; a++) add(scores, scoreError, a, w[row + a]);
                }
                for (int a = 0; a < p; a++) {
                    addProduct(base, error, a, scores[a], scores[a], false);
                    addProduct(base, error, a, scores[a], scoreError[a], true);
                    addProduct(base, error, a, scoreError[a], scoreError[a], false);
                }
            }
            return new double[][] {base, error};
        }
        for (int i = 0; i < n; i++) for (int a = 0; a < p; a++)
            addProduct(base, error, a, w[i * p + a], w[i * p + a], false);
        if (mask.enabled) for (int i = 1; i < n; i++) {
            final int offset = pairIndex(i, 0), limit = offset + i;
            for (int k = mask.bits.nextSetBit(offset); k >= 0 && k < limit; k = mask.bits.nextSetBit(k + 1)) {
                final int j = k - offset;
                for (int a = 0; a < p; a++) addProduct(base, error, a, w[i * p + a], w[j * p + a], true);
            }
        }
        return new double[][] {base, error};
    }

    private static void storePaths(double[] wb, double[] error, double[] w,
            int n, int p, int G, ProtectedPairs mask) {
        final double[][] base = baseline(w, n, p, mask);
        createMatrix("__tmo_paths", p, G);
        for (int a = 0; a < p; a++) {
            for (int g = G - 1; g >= 0; g--) {
                final int bin = a * (G + 1) + g + 1;
                add(base[0], base[1], a, wb[bin]);
                add(base[0], base[1], a, error[bin]);
                final double value = base[0][a] + base[1][a];
                if (!validNumber(value)) throw new ProcedureException("variance path exceeds the finite numeric range");
                Matrix.storeMatrixAt("__tmo_paths", a, g, value);
            }
        }
    }

    private static void storeCounts(long[] bins, int G, long unsupported) {
        createMatrix("__tmo_screened_counts", 1, G);
        long running = 0;
        for (int g = G - 1; g >= 0; g--) {
            running += bins[g + 1];
            Matrix.storeMatrixAt("__tmo_screened_counts", 0, g, running);
        }
        Scalar.setValue("__tmo_unsupported_pairs", (double) unsupported);
        Scalar.setValue("__tmo_counts_api", 1.0);
    }

    private static void checkMemory(int n, int p, int G, int threads,
            boolean endpoint, boolean missing) {
        // Streaming target diagnostics never allocate a common-pair triangle.
        checkMemory(n, p, G, threads, endpoint, missing, false, 0, false, false, 0, false);
    }

    private static void checkMemory(int n, int p, int G, int threads,
            boolean endpoint, boolean missing, boolean masked, long protectedCount, boolean compact) {
        checkMemory(n, p, G, threads, endpoint, missing, masked, protectedCount, compact,
            !missing && !compact && optionalScalar("__tmo_cache_mode", 1.0) != 0.0);
    }

    private static void checkMemory(int n, int p, int G, int threads,
            boolean endpoint, boolean missing, boolean masked, long protectedCount,
            boolean compact, boolean cachePairs) {
        checkMemory(n, p, G, threads, endpoint, missing, masked, protectedCount, compact, cachePairs, 0);
    }

    private static void checkMemory(int n, int p, int G, int threads,
            boolean endpoint, boolean missing, boolean masked, long protectedCount,
            boolean compact, boolean cachePairs, int missingWorkers) {
        checkMemory(n, p, G, threads, endpoint, missing, masked, protectedCount,
            compact, cachePairs, missingWorkers, true);
    }

    private static void checkMemory(int n, int p, int G, int threads,
            boolean endpoint, boolean missing, boolean masked, long protectedCount,
            boolean compact, boolean cachePairs, int missingWorkers, boolean creditAllowed) {
        checkMemory(n, p, G, threads, endpoint, missing, masked, protectedCount,
            compact, cachePairs, missingWorkers, creditAllowed, 0L);
    }

    private static void checkMemory(int n, int p, int G, int threads,
            boolean endpoint, boolean missing, boolean masked, long protectedCount,
            boolean compact, boolean cachePairs, int missingWorkers, boolean creditAllowed,
            long liveInputBytes) {
        final long cells = (long) n * p;
        final long pairs = (long) n * (n - 1) / 2;
        if (cells > Integer.MAX_VALUE - 8L || (long) p * p > Integer.MAX_VALUE - 8L
                || (long) p * (G + 1L) > Integer.MAX_VALUE - 8L
                || (endpoint && pairs - protectedCount > Integer.MAX_VALUE - 8L)
                || (((masked && !compact) || cachePairs) && pairs > Integer.MAX_VALUE - 8L)
                || (masked && !compact && 2 * protectedCount > Integer.MAX_VALUE - 8L)) {
            throw new MemoryLimitException("An exact work array exceeds the Java array limit.");
        }
        final double bytes = (endpoint ? 8.0 * (pairs - protectedCount) : 0.0)
            + (cachePairs ? 8.0 * pairs + 8.0 * cells : 0.0)
            + (compact ? 32.0 * n : masked ? Math.ceil(pairs / 64.0) * 8.0 + 16.0 * protectedCount + 12.0 * n : 0.0)
            + (missing ? 25.0 + 8.0 * missingWorkers : 32.0) * cells
            + 8.0 * p * p + 16.0 * (threads + 1.0) * p * (G + 1.0)
            + 8.0 * n + 64.0 * p + 16.0 * G
            + (optionalScalar("__tmo_need_counts", 0.0) == 1.0 ? 8.0 * (threads + 1.0) * (G + 2.0) : 0.0);
        final Runtime runtime = Runtime.getRuntime();
        final long reserve = Math.max(16L * 1024 * 1024, runtime.maxMemory() / 10);
        long available = runtime.maxMemory() - runtime.totalMemory() + runtime.freeMemory();
        final long reusable = creditAllowed ? reusableBytes(n, p, G, threads, endpoint, missing,
            cachePairs, pairs - protectedCount, missingWorkers) : 0;
        // Credit only exact-size arrays that this phase will actually reuse.
        if (bytes + reserve <= available + reusable + liveInputBytes) return;
        if (workspace.retainedBytes() > 0) {
            workspace.release();
        }
        if (bytes + reserve > available + liveInputBytes) {
            // Only collect under allocation pressure, not on every fit.
            System.gc();
            available = runtime.maxMemory() - runtime.totalMemory() + runtime.freeMemory();
        }
        if (bytes + reserve > available + liveInputBytes) {
            throw new MemoryLimitException("Estimated work arrays require "
                + (long) Math.ceil(bytes / (1024 * 1024)) + " MiB; available heap is "
                + available / (1024 * 1024) + " MiB before the safety reserve.");
        }
    }

    @FunctionalInterface
    private interface Worker { void run(int index); }

    private static void runWorkers(int threads, Worker work) throws Exception {
        if (!reuseWorkers) { runFreshWorkers(threads, work); return; }
        if (threads == 1) { work.run(0); return; }
        if (workerPool == null || workerPoolSize < threads) {
            if (workerPool != null) workerPool.shutdown();
            workerPoolSize = threads;
            workerPool = Executors.newFixedThreadPool(threads, task -> {
                Thread thread = new Thread(task, "tmo-kernel");
                thread.setDaemon(true);
                return thread;
            });
        }
        final CountDownLatch done = new CountDownLatch(threads);
        final Throwable[] failures = new Throwable[threads];
        int submitted = 0;
        try {
            for (int t = 0; t < threads; t++) {
                final int index = t;
                workerPool.execute(() -> {
                    try { work.run(index); }
                    catch (Throwable failure) { failures[index] = failure; }
                    finally { done.countDown(); }
                });
                submitted++;
            }
        } catch (Throwable failure) {
            failures[submitted] = failure;
            // Submission itself can fail after earlier tasks started.
            for (int t = submitted; t < threads; t++) done.countDown();
        }
        // Do not release/reuse arrays while any worker can still access them.
        boolean interrupted = false;
        while (true) {
            try { done.await(); break; }
            catch (InterruptedException failure) { interrupted = true; }
        }
        if (interrupted) Thread.currentThread().interrupt();
        for (Throwable failure : failures) {
            if (failure instanceof Error) throw (Error) failure;
            if (failure instanceof RuntimeException) throw (RuntimeException) failure;
            if (failure != null) throw new Exception(failure);
        }
        if (interrupted) throw new InterruptedException("tmo workers interrupted");
    }

    private static void runFreshWorkers(int threads, Worker work) throws Exception {
        if (threads == 1) { work.run(0); return; }
        final Thread[] pool = new Thread[threads];
        final Throwable[] failures = new Throwable[threads];
        int started = 0;
        boolean interrupted = false;
        try {
            for (int t = 0; t < threads; t++) {
                final int index = t;
                pool[t] = new Thread(() -> {
                    try { work.run(index); }
                    catch (Throwable failure) { failures[index] = failure; }
                });
                pool[t].start();
                started++;
            }
        } finally {
            for (int t = 0; t < started; t++) {
                while (true) {
                    try { pool[t].join(); break; }
                    catch (InterruptedException failure) { interrupted = true; }
                }
            }
            if (interrupted) Thread.currentThread().interrupt();
        }
        for (Throwable failure : failures) {
            if (failure instanceof Error) throw (Error) failure;
            if (failure instanceof RuntimeException) throw (RuntimeException) failure;
            if (failure != null) throw new Exception(failure);
        }
    }

    private static int compute() throws Exception {
        final long preparationStart = System.nanoTime();
        final int cacheMode = (int) optionalScalar("__tmo_cache_mode", 1.0);
        if (cacheMode < 0 || cacheMode > 2) throw new ProcedureException("cache mode must be 0, 1, or 2");
        final boolean reuse = cacheMode != 0;
        reuseWorkers = reuse;
        if (!reuse) workspace.release();
        Scalar.setValue("__tmo_cache_hit", 0.0);
        Scalar.setValue("__tmo_cache_rebuild", 0.0);
        Scalar.setValue("__tmo_cache_fallback_pairs", 0.0);
        final int n = (int) Mata.getMataRowTotal("__tmo_Xs");
        final int p = (int) Mata.getMataColTotal("__tmo_Xs");
        final boolean endpoint = optionalScalar("__tmo_endpoint", 0.0) == 1.0;
        final boolean auxOnly = optionalScalar("__tmo_auxonly", 0.0) == 1.0;
        final double requestedG = endpoint ? Mata.getMataRowTotal("__tmo_etagrid")
                : Mata.getMataRowTotal("__tmo_dgrid");
        if (requestedG != Math.rint(requestedG) || requestedG < 1
                || requestedG > Integer.MAX_VALUE - 1) {
            SFIToolkit.errorln("tmo: grid size must be a positive integer");
            return error(459);
        }
        final int G = (int) requestedG;
        if (n <= 1 || p <= 1 || G < 1) {
            SFIToolkit.errorln("tmo: bad input dimensions");
            return error(459);
        }
        if ((int) Mata.getMataRowTotal("__tmo_W") != n
                || (int) Mata.getMataColTotal("__tmo_W") != p) {
            SFIToolkit.errorln("tmo: __tmo_W must match __tmo_Xs");
            return error(459);
        }

        int requestedThreads = (int) optionalScalar("__tmo_threads", 0.0);
        if (requestedThreads <= 0) requestedThreads = Runtime.getRuntime().availableProcessors();
        final int threads = Math.max(1, Math.min(requestedThreads, n - 1));
        final boolean hasMissing = optionalScalar("__tmo_hasmiss", 0.0) == 1.0;
        final boolean masked = optionalScalar("__tmo_masked", 0.0) == 1.0;
        final boolean compact = masked && Mata.getMataRowTotal("__tmo_Bgroup") > 0;
        final boolean needCounts = optionalScalar("__tmo_need_counts", 0.0) == 1.0;
        long protectedCount = 0;
        if (masked) {
            final long rows = Mata.getMataRowTotal("__tmo_Bpairs");
            if (rows < 0 || Mata.getMataColTotal("__tmo_Bpairs") != 2
                    || rows > (long) n * (n - 1) / 2)
                throw new ProcedureException("__tmo_Bpairs must be a q by 2 strict-lower pair matrix");
            protectedCount = rows;
            if (endpoint && rows == (long) n * (n - 1) / 2) {
                Scalar.setValue("__tmo_calibration_error", 1.0);
                throw new ProcedureException("no eligible off-diagonal pairs remain for IQR calibration");
            }
        }
        // Compact metadata needs O(n) space before the exact eligible count is
        // known. No strict-triangle cache is allocated for this representation.
        checkMemory(n, p, G, threads, endpoint, hasMissing, masked,
            compact ? (long) n * (n - 1) / 2 : protectedCount, compact);
        final ProtectedPairs mask = new ProtectedPairs(n, masked);
        if (compact) checkMemory(n, p, G, threads, endpoint, hasMissing, masked,
            (long) n * (n - 1) / 2 - mask.eligibleCount, true);
        if (endpoint && mask.eligibleCount == 0) {
            Scalar.setValue("__tmo_calibration_error", 1.0);
            throw new ProcedureException("no eligible off-diagonal pairs remain for IQR calibration");
        }

        // Mata stores matrices ROW-major, and getMataReal hands the buffer
        // back in that order, which is already what the inner loop over
        // anchors wants: x[i*p + a].
        final double[] x = Mata.getMataReal("__tmo_Xs");
        final double[] w = Mata.getMataReal("__tmo_W");
        final double[] delta = endpoint ? new double[G] : Mata.getMataReal("__tmo_dgrid");
        final double[] eta = endpoint ? Mata.getMataReal("__tmo_etagrid") : null;

        if (x.length != n*p || w.length != n*p || delta.length != G) {
            SFIToolkit.errorln("tmo: unexpected buffer length");
            return error(459);
        }
        if (endpoint) {
            if ((int) Mata.getMataColTotal("__tmo_etagrid") != 1 || eta.length != G) {
                SFIToolkit.errorln("tmo: __tmo_etagrid must be a G x 1 column vector");
                return error(459);
            }
            for (int g = 0; g < G; g++) {
                if (!validNumber(eta[g]) || !(eta[g] > 0.0) || eta[g] > 1.0
                        || (g > 0 && !(eta[g] > eta[g - 1]))) {
                    SFIToolkit.errorln("tmo: __tmo_etagrid must be strictly increasing within (0,1]");
                    return error(459);
                }
            }
        }

        for (int k = 0; k < x.length; k++) {
            if ((!validNumber(x[k]) && !(hasMissing && Data.isValueMissing(x[k])))
                    || !validNumber(w[k])) {
                SFIToolkit.errorln("tmo: screening residuals and weights must be finite and nonmissing");
                return error(459);
            }
        }
        final boolean corr = optionalScalar("__tmo_corrscale", 0.0) == 1.0;
        if (endpoint && !corr) {
            SFIToolkit.errorln("tmo: endpoint grid requires correlation screens");
            return error(459);
        }
        final int d = p - 1;
        final int first = auxOnly ? 1 : 0;

        if (hasMissing) {
            workspace.invalidate();
            if (!corr) {
                SFIToolkit.errorln("tmo: missing outcomes require correlation screening");
                return error(459);
            }
            return computeMissing(x, w, delta, eta, n, p, G, first, threads, endpoint, mask, needCounts);
        }

        // gram over outcomes, and h_a = gram[.,a]/ss[a]
        final double[] gram = reuse ? (workspace.gram = array(workspace.gram, p * p, true)) : new double[p * p];
        for (int i = 0; i < n; i++) {
            final int base = i * p;
            for (int a = 0; a < p; a++) {
                final double xa = x[base + a];
                if (xa == 0.0) continue;
                for (int b = a; b < p; b++) gram[a * p + b] += xa * x[base + b];
            }
        }
        for (int a = 0; a < p; a++) for (int b = 0; b < a; b++) gram[a * p + b] = gram[b * p + a];
        final double[] ss = reuse ? (workspace.ss = array(workspace.ss, p, false)) : new double[p];
        for (int a = 0; a < p; a++) ss[a] = gram[a * p + a];
        for (int a = 0; a < p; a++) {
            if (!(ss[a] > 0.0) || !validNumber(ss[a])) {
                if (corr) {
                    workspace.invalidate();
                    return computeMissing(x, w, delta, eta, n, p, G, first, threads, endpoint, mask, needCounts);
                }
                SFIToolkit.errorln("tmo: outcome " + (a + 1) + " has no variation");
                return error(459);
            }
        }
        // hh[a] = h_a'h_a
        final double[] hh = reuse ? (workspace.hh = array(workspace.hh, p, false)) : new double[p];
        for (int a = 0; a < p; a++) {
            double acc = 0.0;
            for (int b = first; b < p; b++) {
                final double h = gram[a * p + b] / ss[a];
                acc += h * h;
            }
            hh[a] = acc;
        }
        // z[i,a] = (Xs h_a)_i - x[i,a]*hh[a]/2, row-major
        final double[] z = reuse ? (workspace.z = array(workspace.z, n * p, false)) : new double[n * p];
        for (int i = 0; i < n; i++) {
            final int base = i * p;
            for (int a = 0; a < p; a++) {
                double acc = 0.0;
                final double inv = 1.0 / ss[a];
                for (int b = first; b < p; b++) acc += x[base + b] * (gram[a * p + b] * inv);
                z[base + a] = acc - x[base + a] * hh[a] * 0.5;
            }
        }
        // gd[i] = sum_b x[i,b]^2, so screen_a(i,i) = (gd[i] - 2 z[i,a] x[i,a])/d
        final double[] gd = reuse ? (workspace.gd = array(workspace.gd, n, false)) : new double[n];
        for (int i = 0; i < n; i++) {
            final int base = i * p;
            double acc = 0.0;
            for (int b = first; b < p; b++) acc += x[base + b] * x[base + b];
            gd[i] = acc;
        }
        // per-anchor reciprocal screen sd, for the correlation scale
        final double[] rsd = reuse ? (workspace.rsd = array(workspace.rsd, n * p, false)) : new double[n * p];
        if (corr) {
            for (int i = 0; i < n; i++) {
                final int base = i * p;
                for (int a = 0; a < p; a++) {
                    final double v = (gd[i] - 2.0 * z[base + a] * x[base + a]) / d;
                    if (illConditionedNorm(gd[i], 2.0 * z[base + a] * x[base + a], d)) {
                        workspace.invalidate();
                        return computeMissing(x, w, delta, eta, n, p, G, first, threads, endpoint, mask, needCounts);
                    }
                    rsd[base + a] = 1.0 / Math.sqrt(v);
                }
            }
        }

        final CachePlan cache = prepareCache(x, n, p, first, masked, compact ? 0 : cacheMode, threads);
        if (cache != null && cache.rolling) cache.screenBounds(x, z, rsd, n, p, corr);
        final long calibrationStart = System.nanoTime();
        Scalar.setValue("__tmo_preparation_ms", (calibrationStart - preparationStart) / 1e6);
        if (endpoint) {
            final double logfactor = optionalScalar("__tmo_logfactor", 0.0);
            if (!(logfactor > 0.0)) {
                SFIToolkit.errorln("tmo: endpoint log factor must be finite and positive");
                return error(459);
            }
            final double[] qs = anchorZeroQuantiles(x, z, rsd, gram, ss, n, p, d, first, threads, mask, cache);
            final double sd0 = (qs[1] - qs[0]) / (2.0 * 0.6744897501960817);
            final double v0 = sd0 * sd0;
            final double delta2 = 2.0 * v0 * logfactor;
            if (!(delta2 > 0.0) || !validNumber(delta2)) {
                if (mask.enabled && (!(v0 > 0.0) || !validNumber(v0))) Scalar.setValue("__tmo_calibration_error", 2.0);
                SFIToolkit.errorln("tmo: endpoint grid requires a finite positive IQR variance");
                return error(459);
            }
            for (int g = 0; g < G; g++) delta[g] = Math.sqrt(delta2 * eta[g]);
            Scalar.setValue("__tmo_v0hat", v0);
            Scalar.setValue("__tmo_q25", qs[0]);
            Scalar.setValue("__tmo_q75", qs[1]);
            createMatrix("__tmo_dgrid_out", G, 1);
            for (int g = 0; g < G; g++) Matrix.storeMatrixAt("__tmo_dgrid_out", g, 0, delta[g]);
        } else if (cache != null) {
            final int[][] ranges = splitPairs(n, threads);
            runWorkers(threads, t -> {
                for (int i = ranges[t][0]; i <= ranges[t][1]; i++)
                    for (int j = 0; j < i; j++) cache.fill(x, i, j, p, first);
            });
        }
        for (int g = 0; g < G; g++) {
            if (!(delta[g] > 0.0) || !validNumber(delta[g])
                    || (g > 0 && !(delta[g] > delta[g - 1]))) {
                SFIToolkit.errorln("tmo: thresholds must be finite, positive, and strictly increasing");
                return error(459);
            }
        }

        final long pathsStart = System.nanoTime();
        Scalar.setValue("__tmo_calibration_ms", (pathsStart - calibrationStart) / 1e6);
        // is the grid regular in delta^2? then the bin is one division
        final double d1sq = delta[0] * delta[0];
        double step = G > 1 ? delta[1] * delta[1] - d1sq : 0.0;
        boolean regular = G > 1 && step > 0.0;
        for (int g = 2; g < G && regular; g++) {
            final double expect = d1sq + g * step;
            final double got = delta[g] * delta[g];
            regular = Math.abs(got - expect) <= 1e-10 * (1.0 + Math.abs(expect));
        }

        final boolean isRegular = regular;
        final double sqStep = step;
        final int bins = G + 1;
        if (reuse && (workspace.partial == null || workspace.partial.length != threads))
            workspace.partial = new double[threads][];
        if (reuse && (workspace.partialError == null || workspace.partialError.length != threads))
            workspace.partialError = new double[threads][];
        final double[][] partial = reuse ? workspace.partial : new double[threads][];
        final double[][] partialError = reuse ? workspace.partialError : new double[threads][];
        final long[][] counts = needCounts ? new long[threads][bins] : null;
        final int[][] ranges = splitPairs(n, threads);

        for (int t = 0; t < threads; t++) {
            partial[t] = array(partial[t], p * bins, true);
            partialError[t] = array(partialError[t], p * bins, true);
        }
        if (!corr || delta[0] <= 1.0) runWorkers(threads, t -> accumulate(
            partial[t], partialError[t], needCounts ? counts[t] : null,
            x, w, z, rsd, gram, ss, delta, n, p, d, first, bins, G,
            corr, isRegular, d1sq, sqStep, ranges[t][0], ranges[t][1], mask, cache, t));

        // Fixed-order reduction is reproducible for a fixed thread count.
        final double[] wb = reuse ? (workspace.wb = array(workspace.wb, p * bins, true)) : new double[p * bins];
        final double[] wbError = reuse ? (workspace.wbError = array(workspace.wbError, p * bins, true)) : new double[p * bins];
        for (int t = 0; t < threads; t++) {
            final double[] q = partial[t];
            for (int k = 0; k < wb.length; k++) {
                add(wb, wbError, k, q[k]);
                add(wb, wbError, k, partialError[t][k]);
            }
        }

        storePaths(wb, wbError, w, n, p, G, mask);
        if (needCounts) {
            final long[] totalCounts = new long[bins];
            for (int t = 0; t < threads; t++) for (int g = 0; g < bins; g++)
                totalCounts[g] += counts[t][g];
            storeCounts(totalCounts, G, 0);
        }
        Scalar.setValue("__tmo_regular", regular ? 1 : 0);
        Scalar.setValue("__tmo_threads_used", threads);
        if (cache != null) cache.commit(x, n, p);
        Scalar.setValue("__tmo_paths_ms", (System.nanoTime() - pathsStart) / 1e6);
        return 0;
    }

    /** Pair sums, with an opt-in, input-verified rolling simulation context. */
    private static final class CachePlan {
        final double[] common, previousX;
        final boolean rolling, eligible;
        final double context, window, maximumX, error;
        final int shifts;
        final long[] fallbacks;
        double[] maxX, maxZ, maxRsd;

        CachePlan(double[] common, double[] previousX, boolean rolling,
                boolean eligible, double context, double window, double maximumX,
                int shifts, double error, int threads) {
            this.common = common; this.previousX = previousX;
            this.rolling = rolling; this.eligible = eligible;
            this.context = context; this.window = window; this.maximumX = maximumX;
            this.shifts = shifts; this.error = error;
            fallbacks = new long[threads];
        }

        void screenBounds(double[] x, double[] z, double[] rsd, int n, int p, boolean corr) {
            maxX = new double[n]; maxZ = new double[n]; maxRsd = new double[n];
            for (int i = 0; i < n; i++) for (int a = 0; a < p; a++) {
                maxX[i] = Math.max(maxX[i], Math.abs(x[i * p + a]));
                maxZ[i] = Math.max(maxZ[i], Math.abs(z[i * p + a]));
                maxRsd[i] = corr ? Math.max(maxRsd[i], rsd[i * p + a]) : 1.0;
            }
        }

        double scoreError(int i, int j, double value, int d) {
            return screenError(error, value, Math.nextUp(maxZ[i] * maxX[j]),
                Math.nextUp(maxX[i] * maxZ[j]), 1.0 / d, maxRsd[i], maxRsd[j]);
        }

        double fill(double[] x, int i, int j, int p, int first) {
            final int k = pairIndex(i, j), bi = i * p, bj = j * p;
            if (rolling) common[k] = (common[k] - previousX[bi] * previousX[bj])
                + x[bi + p - 1] * x[bj + p - 1];
            else common[k] = dot(x, bi, bj, p, first);
            return common[k];
        }

        void commit(double[] x, int n, int p) {
            workspace.rollingValid = eligible;
            workspace.previousX = eligible ? x : null;
            workspace.previousN = n; workspace.previousP = p;
            workspace.context = context; workspace.window = window;
            workspace.maximumX = maximumX; workspace.shifts = shifts;
            if (rolling) workspace.hits++; else workspace.rebuilds++;
            long count = 0;
            for (long value : fallbacks) count += value;
            Scalar.setValue("__tmo_cache_hit", rolling ? 1 : 0);
            Scalar.setValue("__tmo_cache_rebuild", rolling ? 0 : 1);
            Scalar.setValue("__tmo_cache_hits_total", workspace.hits);
            Scalar.setValue("__tmo_cache_rebuilds_total", workspace.rebuilds);
            Scalar.setValue("__tmo_cache_fallback_pairs", count);
        }
    }

    private static CachePlan prepareCache(double[] x, int n, int p, int first,
            boolean masked, int mode, int threads) {
        if (mode == 0) { workspace.invalidate(); return null; }
        final long pairs = (long) n * (n - 1) / 2;
        if (pairs > Integer.MAX_VALUE - 8L) throw new MemoryLimitException("pair cache exceeds the Java array limit");
        final double context = optionalScalar("__tmo_sim_context", Double.NaN);
        final double window = optionalScalar("__tmo_sim_window", Double.NaN);
        final boolean eligible = mode == 2 && first == 0 && !masked && p <= 4096
            && optionalScalar("__tmo_sim_cache", 0.0) == 1.0
            && Double.isFinite(context) && context == Math.rint(context)
            && Double.isFinite(window) && window >= 1.0 && window < 0x1.0p53
            && window == Math.rint(window);
        boolean rolling = eligible && workspace.rollingValid
            && workspace.previousN == n && workspace.previousP == p
            && workspace.context == context && workspace.window + 1.0 == window
            && workspace.shifts < 32 && workspace.common != null
            && workspace.common.length == (int) pairs;
        double maximum = 0.0;
        for (double value : x) maximum = Math.max(maximum, Math.abs(value));
        if (rolling) {
            // Verifies source values, column order, sample order and every
            // overlapping cell; IDs alone never authorize reuse.
            for (int i = 0; i < n && rolling; i++) for (int a = 0; a < p - 1; a++) {
                if (Double.doubleToRawLongBits(x[i * p + a])
                        != Double.doubleToRawLongBits(workspace.previousX[i * p + a + 1])) {
                    rolling = false;
                    break;
                }
            }
        }
        final int shifts = rolling ? workspace.shifts + 1 : 0;
        if (rolling) maximum = Math.max(maximum, workspace.maximumX);
        final double length = p + 2.0 * shifts;
        // A fresh p-term dot and a chain of p+2*shifts rounded products/adds
        // differ by at most this conservative forward-error bound. For
        // p<=4096, shifts<=32, 16*ulp(1)*length^2 dominates the two standard
        // gamma_k bounds, product rounding, and bound-evaluation rounding.
        // The additive term covers gradual underflow. Nonfinite bounds
        // disable rolling. No compensated/reordered dot defines the reference.
        final double productBound = Math.nextUp(maximum * maximum);
        final double error = Math.nextUp(16.0 * Math.ulp(1.0) * length * length * productBound
            + 16.0 * length * Double.MIN_VALUE);
        if (!Double.isFinite(error)) rolling = false;
        final double[] oldX = rolling ? workspace.previousX : null;
        workspace.rollingValid = false; // failures cannot commit a partly updated triangle
        workspace.common = array(workspace.common, (int) pairs, false);
        return new CachePlan(workspace.common, oldX, rolling, eligible,
            context, window, maximum, rolling ? shifts : 0, rolling ? error : 0.0, threads);
    }

    private static double dot(double[] x, int bi, int bj, int p, int first) {
        double common = 0.0;
        for (int b = first; b < p; b++) common += x[bi + b] * x[bj + b];
        return common;
    }

    /** Conservative difference from evaluating the same screen with a fresh dot. */
    private static double screenError(double commonError, double common,
            double left, double right, double invd, double rsi, double rsj) {
        double error = commonError + 32.0 * Math.ulp(1.0)
            * (Math.abs(common) + Math.abs(left) + Math.abs(right)) + 32.0 * Double.MIN_VALUE;
        error = Math.nextUp(error);
        error = Math.nextUp(error * Math.abs(invd) + 32.0 * Double.MIN_VALUE);
        error = Math.nextUp(error * Math.abs(rsi) + 32.0 * Double.MIN_VALUE);
        return Math.nextUp(error * Math.abs(rsj) + 32.0 * Double.MIN_VALUE);
    }

    private static boolean uncertainBin(double score, double error, int bin, double[] delta) {
        final double value = Math.abs(score);
        return !Double.isFinite(error)
            || (bin > 0 && value - error <= delta[bin - 1])
            || (bin < delta.length && value + error >= delta[bin]);
    }

    /** Split the i>j pair triangle into contiguous i-blocks of even work. */
    private static int[][] splitPairs(int n, int threads) {
        final long total = (long) n * (n - 1) / 2;
        final int[][] out = new int[threads][2];
        int prev = 1;
        for (int t = 0; t < threads; t++) {
            final long target = total * (t + 1) / threads;
            int i = prev;
            long acc = (long) (prev - 1) * prev / 2;   // pairs with row < prev
            // i is an EXCLUSIVE boundary. Empty ranges are allowed once an
            // earlier block has already crossed this worker's pair target.
            while (i < n && acc < target) { acc += i; i++; }
            out[t][0] = prev;
            out[t][1] = (t == threads - 1) ? n - 1 : i - 1;
            prev = out[t][1] + 1;
        }
        out[threads - 1][1] = n - 1;
        return out;
    }

    private static void accumulate(double[] wb, double[] error, long[] counts,
            double[] x, double[] w, double[] z,
            double[] rsd, double[] gram, double[] ss, double[] delta, int n, int p, int d, int first,
            int bins, int G, boolean corr, boolean regular, double d1sq,
            double step, int iLo, int iHi, ProtectedPairs mask, CachePlan cache, int worker) {
        if (cache != null && !cache.rolling && !mask.enabled && corr && regular
                && Double.isFinite(1.0 / step)) {
            accumulateCompleteRegular(wb, error, counts, x, w, z, rsd, gram, ss, delta,
                p, d, first, bins, G, d1sq, 1.0 / step, iLo, iHi, cache.common);
            return;
        }
        final double invd = 1.0 / d;
        for (int i = iLo; i <= iHi; i++) {
            final int bi = i * p;
            for (int j = 0; j < i; j++) {
                if (mask.contains(i, j)) continue;
                final int bj = j * p;
                double common = cache == null ? dot(x, bi, bj, p, first) : cache.common[pairIndex(i, j)];
                boolean approximate = cache != null && cache.rolling;
                final double uncertainty = approximate ? cache.scoreError(i, j, common, d) : 0.0;
                for (int a = 0; a < p; a++) {
                    if ((a != 0 || counts == null) && (w[bi + a] == 0.0 || w[bj + a] == 0.0)) continue;
                    final double left = z[bi + a] * x[bj + a], right = x[bi + a] * z[bj + a];
                    double s = (common - left - right) * invd;
                    if (corr) s = boundCorrelation(s * rsd[bi + a] * rsd[bj + a]);
                    int bin = binOf(s, delta, G, regular, d1sq, step);
                    if (approximate && uncertainBin(s, uncertainty, bin, delta)) {
                        common = dot(x, bi, bj, p, first);
                        cache.common[pairIndex(i, j)] = common;
                        cache.fallbacks[worker]++;
                        approximate = false;
                        s = (common - left - right) * invd;
                        if (corr) s = boundCorrelation(s * rsd[bi + a] * rsd[bj + a]);
                        bin = binOf(s, delta, G, regular, d1sq, step);
                    }
                    if (corr && (illConditionedPair(common, left, right) || nearPerfect(s))) {
                        s = explicitCompletePair(x, gram, ss, bi, bj, p, a, first);
                        bin = binOf(s, delta, G, regular, d1sq, step);
                    }
                    if (a == 0 && counts != null) counts[bin]++;
                    if (bin > 0) addProduct(wb, error, a * bins + bin, w[bi + a], w[bj + a], true);
                }
            }
        }
    }

    /** Complete-data, unmasked correlation path with literal threshold comparisons. */
    private static void accumulateCompleteRegular(double[] wb, double[] error, long[] counts,
            double[] x, double[] w, double[] z, double[] rsd, double[] gram, double[] ss, double[] delta,
            int p, int d, int first, int bins, int G, double firstSquared, double inverseStep,
            int iLo, int iHi, double[] commonPairs) {
        final double invd = 1.0 / d, minimum = delta[0], maximum = delta[G - 1];
        final double[] scores = new double[p];
        for (int i = iLo; i <= iHi; i++) {
            final int bi = i * p, row = pairIndex(i, 0);
            for (int j = 0; j < i; j++) {
                final int bj = j * p;
                final double common = commonPairs[row + j];
                for (int a = 0; a < p; a++) {
                    if ((a != 0 || counts == null) && (w[bi + a] == 0.0 || w[bj + a] == 0.0)) continue;
                    final double left = z[bi + a] * x[bj + a], right = x[bi + a] * z[bj + a];
                    final double fast = ((common - left - right) * invd) * rsd[bi + a] * rsd[bj + a];
                    final double s = illConditionedPair(common, left, right) || nearPerfect(fast)
                        ? explicitCompletePair(x, gram, ss, bi, bj, p, a, first)
                        : boundCorrelation(fast);
                    scores[a] = s;
                }
                for (int a = 0; a < p; a++) {
                    if ((a != 0 || counts == null) && (w[bi + a] == 0.0 || w[bj + a] == 0.0)) continue;
                    final double s = scores[a], absolute = Math.abs(s);
                    int bin;
                    if (absolute < minimum) bin = 0;
                    else if (absolute >= maximum) bin = G;
                    else {
                        // Multiplication only proposes a bin; literal cutoffs decide it.
                        bin = Math.min(G, Math.max(0, 1 + (int) Math.floor((s * s - firstSquared) * inverseStep)));
                        while (bin > 0 && delta[bin - 1] > absolute) bin--;
                        while (bin < G && delta[bin] <= absolute) bin++;
                    }
                    if (a == 0 && counts != null) counts[bin]++;
                    if (bin > 0) addProduct(wb, error, a * bins + bin, w[bi + a], w[bj + a], true);
                }
            }
        }
    }

    /** Number of grid thresholds that |score| clears. */
    private static int binOf(double score, double[] delta, int G,
            boolean regular, double d1sq, double step) {
        final double a = Math.abs(score);
        if (regular) {
            if (a < delta[0]) return 0;
            if (a >= delta[G - 1]) return G;
            final double sq = score * score;
            int cleared = Math.min(G, Math.max(0, 1 + (int) Math.floor((sq - d1sq) / step)));
            // The regular-grid calculation is only a fast candidate. Compare
            // with actual thresholds to preserve the literal >= boundary,
            // including one-ulp neighbors and nearly regular user grids.
            while (cleared > 0 && delta[cleared - 1] > a) cleared--;
            while (cleared < G && delta[cleared] <= a) cleared++;
            return cleared;
        }
        int lo = 0, hi = G;
        while (lo < hi) {
            final int mid = lo + (hi - lo) / 2;
            if (delta[mid] <= a) lo = mid + 1; else hi = mid;
        }
        return lo;
    }

    private static boolean validNumber(double x) {
        return Double.isFinite(x) && !Data.isValueMissing(x);
    }

    // Enforce the correlation range, without a threshold-comparison tolerance.
    private static double boundCorrelation(double score) {
        return Math.max(-1.0, Math.min(1.0, score));
    }

    // This gate selects stable arithmetic; it never declares a pair perfect.
    private static boolean nearPerfect(double score) {
        return Math.abs(score) >= 1.0 - Math.sqrt(Math.ulp(1.0));
    }

    // This selects an arithmetic evaluation, never a threshold tie tolerance.
    private static boolean illConditionedNorm(double total, double subtract, int d) {
        final double numerator = total - subtract;
        return !validNumber(numerator) || numerator <= Math.sqrt(Math.ulp(1.0))
            * (Math.abs(total) + Math.abs(subtract)) || !(numerator / d > 0.0);
    }

    private static boolean illConditionedPair(double common, double left, double right) {
        return !Double.isFinite(common - left - right) || Math.abs(common - left - right)
            <= Math.sqrt(Math.ulp(1.0)) * (Math.abs(common) + Math.abs(left) + Math.abs(right));
    }

    private static double explicitCompletePair(double[] x, double[] gram, double[] ss,
            int bi, int bj, int p, int anchor, int first) {
        double maxI = 0.0, maxJ = 0.0, singleI = 0.0, singleJ = 0.0;
        int shared = 0;
        for (int b = first; b < p; b++) {
            if (b == anchor) continue;
            final double h = gram[anchor * p + b] / ss[anchor];
            singleI = x[bi + b] - x[bi + anchor] * h;
            singleJ = x[bj + b] - x[bj + anchor] * h;
            maxI = Math.max(maxI, Math.abs(singleI)); maxJ = Math.max(maxJ, Math.abs(singleJ)); shared++;
        }
        if (!(maxI > 0.0) || !(maxJ > 0.0)) return 0.0;
        if (shared == 1) return Math.signum(singleI) * Math.signum(singleJ);
        double num = 0.0, si = 0.0, sj = 0.0;
        boolean same = true, opposite = true;
        for (int b = first; b < p; b++) {
            if (b == anchor) continue;
            final double h = gram[anchor * p + b] / ss[anchor];
            final double ei = (x[bi + b] - x[bi + anchor] * h) / maxI;
            final double ej = (x[bj + b] - x[bj + anchor] * h) / maxJ;
            same &= ei == ej; opposite &= ei == -ej;
            num += ei * ej; si += ei * ei; sj += ej * ej;
        }
        final double score = num / Math.sqrt(si) / Math.sqrt(sj);
        if (same) return 1.0;
        if (opposite) return -1.0;
        return validNumber(score) ? boundCorrelation(score) : 0.0;
    }

    private static double optionalScalar(String name, double fallback) {
        final Double value = Scalar.getValue(name);
        return value == null || !validNumber(value.doubleValue()) ? fallback : value.doubleValue();
    }

    private static int computeMissing(double[] x, double[] w, double[] delta,
            double[] eta, int n, int p, int G, int first, int threads,
            boolean endpoint, ProtectedPairs mask, boolean needCounts) throws Exception {
        Scalar.setValue("__tmo_explicit_screen", 1.0);
        final long preparationStart = activeRunStart;
        checkMemory(n, p, G, threads, endpoint, true, mask.enabled,
            (long) n * (n - 1) / 2 - mask.eligibleCount, mask.groups != null, false, 0, true, 16L * x.length);
        if (reuseWorkers && (workspace.observed == null || workspace.observed.length != x.length))
            workspace.observed = new boolean[x.length];
        final boolean[] observed = reuseWorkers ? workspace.observed : new boolean[x.length];
        for (int k = 0; k < x.length; k++) {
            observed[k] = validNumber(x[k]);
            if (!observed[k]) x[k] = 0.0;
        }
        final long[] unsupported = needCounts ? new long[threads] : null;
        final long calibrationStart = System.nanoTime();
        Scalar.setValue("__tmo_preparation_ms", (calibrationStart - preparationStart) / 1e6);
        if (endpoint || needCounts) {
            final double[] residual = reuseWorkers
                ? (workspace.missingResidual = array(workspace.missingResidual, x.length, false)) : new double[x.length];
            final boolean supported = missingResiduals(x, observed, residual, n, p, 0, first);
            final double[] lower = endpoint ? (reuseWorkers
                ? (workspace.lower = array(workspace.lower, (int) mask.eligibleCount, true)) : new double[(int) mask.eligibleCount]) : null;
            final int[][] ranges = splitPairs(n, threads);
            if (!supported && needCounts) unsupported[0] = (long) n * (n - 1) / 2;
            if (supported) runWorkers(threads, t -> {
                final int[] details = needCounts ? new int[2] : null;
                for (int i = ranges[t][0]; i <= ranges[t][1]; i++) {
                    int offset = mask.offset(i);
                    for (int j = 0; j < i; j++) {
                        final boolean protectedPair = mask.contains(i, j);
                        if (protectedPair && !needCounts) continue;
                        final double score = missingScore(residual, observed, i, j, p, 0, first, details);
                        if (needCounts && details[1] == 0) unsupported[t]++;
                        if (endpoint && !protectedPair) lower[offset++] = score;
                    }
                }
            });
            if (endpoint) {
            final double[] qs = calibrationQuartiles(lower, n, mask.enabled);
            final double sd0 = (qs[1] - qs[0]) / (2.0 * 0.6744897501960817);
            final double v0 = sd0 * sd0;
            final double logfactor = optionalScalar("__tmo_logfactor", 0.0);
            final double delta2 = 2.0 * v0 * logfactor;
            if (!(delta2 > 0.0) || !validNumber(delta2)) {
                if (mask.enabled && (!(v0 > 0.0) || !validNumber(v0))) Scalar.setValue("__tmo_calibration_error", 2.0);
                SFIToolkit.errorln("tmo: endpoint grid requires a finite positive IQR variance and log factor");
                return error(459);
            }
            for (int g = 0; g < G; g++) delta[g] = Math.sqrt(delta2 * eta[g]);
            Scalar.setValue("__tmo_v0hat", v0);
            Scalar.setValue("__tmo_q25", qs[0]);
            Scalar.setValue("__tmo_q75", qs[1]);
            createMatrix("__tmo_dgrid_out", G, 1);
            for (int g = 0; g < G; g++) Matrix.storeMatrixAt("__tmo_dgrid_out", g, 0, delta[g]);
            }
        }
        for (int g = 0; g < G; g++) {
            if (!(delta[g] > 0.0) || !validNumber(delta[g])
                    || (g > 0 && !(delta[g] > delta[g - 1]))) {
                SFIToolkit.errorln("tmo: thresholds must be finite, positive, and strictly increasing");
                return error(459);
            }
        }
        final double d1sq = delta[0] * delta[0];
        final long pathsStart = System.nanoTime();
        Scalar.setValue("__tmo_calibration_ms", (pathsStart - calibrationStart) / 1e6);
        final double step = G > 1 ? delta[1] * delta[1] - d1sq : 0.0;
        boolean regular = G > 1 && step > 0.0;
        for (int g = 2; g < G && regular; g++) {
            final double expected = d1sq + g * step;
            regular = Math.abs(delta[g] * delta[g] - expected)
                <= 1e-10 * (1.0 + Math.abs(expected));
        }
        final boolean isRegular = regular;
        final int bins = G + 1;
        int fittedWorkers = Math.min(threads, p);
        if (delta[0] <= 1.0) {
            while (true) {
                try {
                    checkMemory(n, p, G, threads, endpoint, true, mask.enabled,
                        (long) n * (n - 1) / 2 - mask.eligibleCount, mask.groups != null, false, fittedWorkers, true, 16L * x.length);
                    break;
                } catch (MemoryLimitException insufficient) {
                    if (fittedWorkers == 1) throw insufficient;
                    fittedWorkers = Math.max(1, fittedWorkers / 2);
                }
            }
        }
        final double[] wb = reuseWorkers ? (workspace.wb = array(workspace.wb, p * bins, true)) : new double[p * bins];
        final double[] wbError = reuseWorkers ? (workspace.wbError = array(workspace.wbError, p * bins, true)) : new double[p * bins];
        final long[] counts = needCounts ? new long[bins] : null;
        final int workers = fittedWorkers;
        if (reuseWorkers && delta[0] <= 1.0) {
            if (workspace.missingWorkers == null || workspace.missingWorkers.length != workers)
                workspace.missingWorkers = new double[workers][];
            for (int t = 0; t < workers; t++)
                workspace.missingWorkers[t] = array(workspace.missingWorkers[t], x.length, false);
        }
        if (delta[0] <= 1.0) runWorkers(workers, t -> {
            final double[] residual = reuseWorkers ? workspace.missingWorkers[t] : new double[x.length];
            for (int a = t; a < p; a += workers) {
                if (!missingResiduals(x, observed, residual, n, p, a, first)) continue;
                final int offset = a * bins;
                for (int i = 1; i < n; i++) {
                    for (int j = 0; j < i; j++) {
                        if (mask.contains(i, j)) continue;
                        if ((a != 0 || !needCounts) && (w[i * p + a] == 0.0 || w[j * p + a] == 0.0)) continue;
                        final double score = missingScore(residual, observed, i, j, p, a, first);
                        final int bin = binOf(score, delta, G, isRegular, d1sq, step);
                        if (a == 0 && needCounts) counts[bin]++;
                        if (bin > 0) addProduct(wb, wbError, offset + bin, w[i * p + a], w[j * p + a], true);
                    }
                }
            }
        });
        storePaths(wb, wbError, w, n, p, G, mask);
        if (needCounts) {
            long totalUnsupported = 0;
            for (long value : unsupported) totalUnsupported += value;
            storeCounts(counts, G, totalUnsupported);
        }
        Scalar.setValue("__tmo_regular", regular ? 1 : 0);
        Scalar.setValue("__tmo_threads_used", workers);
        Scalar.setValue("__tmo_paths_ms", (System.nanoTime() - pathsStart) / 1e6);
        return 0;
    }

    // A slope is fitted only on observations shared by its outcome and anchor.
    // An undefined projection removes that auxiliary column from this anchor's
    // support. NaN marks omitted columns; raw observation masks stay immutable.
    private static boolean missingResiduals(double[] x, boolean[] observed,
            double[] residual, int n, int p, int anchor, int first) {
        boolean any = false;
        for (int b = 0; b < p; b++) {
            if (b < first || b == anchor) {
                for (int i = 0; i < n; i++) residual[i * p + b] = 0.0;
                continue;
            }
            double numerator = 0.0, denominator = 0.0, maxA = 0.0, maxB = 0.0;
            for (int i = 0; i < n; i++) {
                final int base = i * p;
                if (observed[base + anchor] && observed[base + b]) {
                    final double a = x[base + anchor], y = x[base + b];
                    numerator += a * y;
                    denominator += a * a;
                    maxA = Math.max(maxA, Math.abs(a));
                    maxB = Math.max(maxB, Math.abs(y));
                }
            }
            if (!(maxA > 0.0)) {
                for (int i = 0; i < n; i++) residual[i * p + b] = Double.NaN;
                continue;
            }
            any = true;
            double slope = numerator / denominator;
            boolean scaled = !validNumber(slope) || !validNumber(numerator)
                || !validNumber(denominator) || denominator < Double.MIN_NORMAL
                || (maxB > 0.0 && (!(maxA * maxB >= Double.MIN_NORMAL)
                    || !Double.isFinite(maxA * maxB)))
                || (numerator != 0.0 && Math.abs(numerator) < Double.MIN_NORMAL);
            if (!scaled) {
                for (int i = 0; i < n; i++) {
                    final int base = i * p;
                    final double e = observed[base + anchor] && observed[base + b]
                        ? x[base + b] - x[base + anchor] * slope : 0.0;
                    residual[base + b] = e;
                    if (!validNumber(e)) scaled = true;
                }
            }
            if (scaled) {
                // Maxima are on the anchor/outcome overlap, never unshared rows.
                double num = 0.0, den = 0.0;
                for (int i = 0; i < n; i++) {
                    final int base = i * p;
                    if (!observed[base + anchor] || !observed[base + b]) continue;
                    final double a = x[base + anchor] / maxA;
                    final double y = maxB > 0.0 ? x[base + b] / maxB : 0.0;
                    final double product = a * y;
                    if ((x[base + anchor] != 0.0 && Math.abs(a) < Double.MIN_NORMAL)
                            || (x[base + b] != 0.0 && Math.abs(y) < Double.MIN_NORMAL)
                            || (a != 0.0 && y != 0.0 && Math.abs(product) < Double.MIN_NORMAL))
                        throw new ProcedureException("auxiliary projection loses precision from extreme within-outcome "
                            + "dynamic range; inspect the observations or transform the affected outcome");
                    num += product;
                    den += a * a;
                }
                slope = num / den;
                for (int i = 0; i < n; i++) {
                    final int base = i * p;
                    final double e = observed[base + anchor] && observed[base + b] && maxB > 0.0
                        ? (x[base + b] / maxB - (x[base + anchor] / maxA) * slope) * maxB : 0.0;
                    if (!validNumber(e)) throw new ProcedureException("auxiliary projection overflow; rescale outcome units");
                    residual[base + b] = e;
                }
            }
        }
        return any;
    }

    private static double missingScore(double[] residual, boolean[] observed,
            int i, int j, int p, int anchor, int first) {
        return missingScore(residual, observed, i, j, p, anchor, first, null);
    }

    private static double missingScore(double[] residual, boolean[] observed,
            int i, int j, int p, int anchor, int first, int[] details) {
        final int bi = i * p, bj = j * p;
        if (details != null) { details[0] = 0; details[1] = 0; }
        if (!observed[bi + anchor] || !observed[bj + anchor]) return 0.0;
        double maxI = 0.0, maxJ = 0.0;
        double numerator = 0.0, ssI = 0.0, ssJ = 0.0;
        double singleI = 0.0, singleJ = 0.0;
        int shared = 0;
        for (int b = first; b < p; b++) {
            if (b == anchor || !observed[bi + b] || !observed[bj + b]) continue;
            final double ei = residual[bi + b], ej = residual[bj + b];
            if (!Double.isFinite(ei) || !Double.isFinite(ej)) continue;
            shared++;
            singleI = ei;
            singleJ = ej;
            maxI = Math.max(maxI, Math.abs(ei));
            maxJ = Math.max(maxJ, Math.abs(ej));
            numerator += ei * ej;
            ssI += ei * ei;
            ssJ += ej * ej;
        }
        if (details != null) details[0] = shared;
        if (!(maxI > 0.0) || !(maxJ > 0.0) || !validNumber(maxI) || !validNumber(maxJ)) return 0.0;
        if (details != null) details[1] = 1;
        if (shared == 1) return Math.signum(singleI) * Math.signum(singleJ);
        // Separate square roots avoid overflow in ssI*ssJ on ordinary data.
        // Exceptional product/square ranges use shared-support scaling below.
        if (validNumber(numerator) && validNumber(ssI) && validNumber(ssJ)
                && ssI >= Double.MIN_NORMAL && ssJ >= Double.MIN_NORMAL) {
            final double partial = numerator / Math.sqrt(ssI);
            final double score = partial / Math.sqrt(ssJ);
            if (validNumber(score) && (numerator == 0.0 || (partial != 0.0 && score != 0.0))
                    && (numerator == 0.0 || Math.abs(numerator) >= Double.MIN_NORMAL)
                    && !nearPerfect(score))
                return boundCorrelation(score);
        }
        // Scale on this pair's shared support, not the whole row: unshared
        // huge entries must not underflow its shared small residuals.
        numerator = ssI = ssJ = 0.0;
        boolean same = true, opposite = true;
        for (int b = first; b < p; b++) {
            if (b == anchor || !observed[bi + b] || !observed[bj + b]) continue;
            if (!Double.isFinite(residual[bi + b]) || !Double.isFinite(residual[bj + b])) continue;
            final double ei = residual[bi + b] / maxI, ej = residual[bj + b] / maxJ;
            same &= ei == ej; opposite &= ei == -ej;
            numerator += ei * ej;
            ssI += ei * ei;
            ssJ += ej * ej;
        }
        final double score = numerator / Math.sqrt(ssI) / Math.sqrt(ssJ);
        if (same) return 1.0;
        if (opposite) return -1.0;
        return validNumber(score) ? boundCorrelation(score) : 0.0;
    }

    /** Exact inverse-ECDF quartiles under the requested calibration population. */
    private static double[] anchorZeroQuantiles(double[] x, double[] z,
            double[] rsd, double[] gram, double[] ss, int n, int p, int d, int first, int threads,
            ProtectedPairs mask, CachePlan cache) throws Exception {
        final long pairs = mask.eligibleCount;
        if (pairs > Integer.MAX_VALUE - 8) throw new IllegalArgumentException("too many pairs for IQR");
        final double[] lower = reuseWorkers
            ? (workspace.lower = array(workspace.lower, (int) pairs, false)) : new double[(int) pairs];
        final double[] errors = new double[threads];
        final long[][] rankCounts = reuseWorkers && !mask.enabled ? new long[threads][2] : null;
        final int[][] ranges = splitPairs(n, threads);
        runWorkers(threads, t -> {
            long belowOne = 0, equalOne = 0;
            for (int i = ranges[t][0]; i <= ranges[t][1]; i++) {
                final int bi = i * p;
                int offset = mask.offset(i);
                for (int j = 0; j < i; j++) {
                    if (mask.contains(i, j)) continue;
                    final int bj = j * p;
                    final double common = cache == null ? dot(x, bi, bj, p, first)
                        : cache.fill(x, i, j, p, first);
                    final double left = z[bi] * x[bj], right = x[bi] * z[bj];
                    final double fast = ((common - left - right) * (1.0 / d)) * rsd[bi] * rsd[bj];
                    final double s = illConditionedPair(common, left, right) || nearPerfect(fast)
                        ? explicitCompletePair(x, gram, ss, bi, bj, p, 0, first)
                        : fast;
                    if (!validNumber(s)) throw new IllegalArgumentException("nonfinite anchor-0 correlation");
                    final double bounded = boundCorrelation(s);
                    lower[offset++] = bounded;
                    if (rankCounts != null) {
                        if (bounded < 1.0) belowOne++; else if (bounded == 1.0) equalOne++;
                    }
                    if (cache != null && cache.rolling) errors[t] = Math.max(errors[t],
                        cache.scoreError(i, j, common, d));
                }
            }
            if (rankCounts != null) { rankCounts[t][0] = belowOne; rankCounts[t][1] = equalOne; }
        });
        final double[] preliminary = calibrationQuartiles(lower, n, mask.enabled, rankCounts);
        if (cache == null || !cache.rolling) return preliminary;
        double maximumError = 0.0;
        for (double value : errors) maximumError = Math.max(maximumError, value);
        final double band = Math.nextUp(2.0 * maximumError);
        // If each score moves by <=E, its kth order statistic moves by <=E.
        // Thus every possible exact kth candidate lies within 2E of the
        // approximate kth value. Recompute that entire band, including ties.
        // Outside scores remain strictly on the same side of either rank.
        // This second pass also restores pair order after in-place selection.
        runWorkers(threads, t -> {
            long belowOne = 0, equalOne = 0;
            for (int i = ranges[t][0]; i <= ranges[t][1]; i++) {
                final int bi = i * p;
                int offset = mask.offset(i);
                for (int j = 0; j < i; j++) {
                    if (mask.contains(i, j)) continue;
                    final int bj = j * p, k = pairIndex(i, j);
                    double common = cache.common[k];
                    double score = boundCorrelation(((common - z[bi] * x[bj] - x[bi] * z[bj])
                        * (1.0 / d)) * rsd[bi] * rsd[bj]);
                    if (illConditionedPair(common, z[bi] * x[bj], x[bi] * z[bj]) || nearPerfect(score))
                        score = explicitCompletePair(x, gram, ss, bi, bj, p, 0, first);
                    if (!Double.isFinite(band) || Math.abs(score - preliminary[0]) <= band
                            || Math.abs(score - preliminary[1]) <= band || Math.abs(score - 1.0) <= band) {
                        common = dot(x, bi, bj, p, first);
                        cache.common[k] = common;
                        cache.fallbacks[t]++;
                        score = boundCorrelation(((common - z[bi] * x[bj] - x[bi] * z[bj])
                            * (1.0 / d)) * rsd[bi] * rsd[bj]);
                        if (illConditionedPair(common, z[bi] * x[bj], x[bi] * z[bj]) || nearPerfect(score))
                            score = explicitCompletePair(x, gram, ss, bi, bj, p, 0, first);
                    }
                    lower[offset++] = score;
                    if (rankCounts != null) {
                        if (score < 1.0) belowOne++; else if (score == 1.0) equalOne++;
                    }
                }
            }
            if (rankCounts != null) { rankCounts[t][0] = belowOne; rankCounts[t][1] = equalOne; }
        });
        return calibrationQuartiles(lower, n, mask.enabled, rankCounts);
    }

    private static double[] calibrationQuartiles(double[] lower, int n,
            boolean eligibleOnly, long[][] rankCounts) {
        if (rankCounts == null) return calibrationQuartiles(lower, n, eligibleOnly);
        long belowOne = 0, equalOne = 0;
        for (long[] counts : rankCounts) { belowOne += counts[0]; equalOne += counts[1]; }
        return allPairQuartiles(lower, n, belowOne, equalOne);
    }

    private static double[] calibrationQuartiles(double[] lower, int n, boolean eligibleOnly) {
        if (!eligibleOnly) return allPairQuartiles(lower, n);
        if (lower.length == 0) {
            Scalar.setValue("__tmo_calibration_error", 1.0);
            throw new ProcedureException("no eligible off-diagonal pairs remain for IQR calibration");
        }
        final long count = 2L * lower.length;
        final long rank25 = (count + 3) / 4;
        final long rank75 = (3 * count + 3) / 4;
        return selectTwo(lower, (int) ((rank25 - 1) / 2), (int) ((rank75 - 1) / 2));
    }

    private static double[] allPairQuartiles(double[] lower, int n) {
        long belowOne = 0, equalOne = 0;
        for (double value : lower) {
            if (value < 1.0) belowOne++;
            else if (value == 1.0) equalOne++;
        }
        return allPairQuartiles(lower, n, belowOne, equalOne);
    }

    private static double[] allPairQuartiles(double[] lower, int n, long belowOne, long equalOne) {
        final long count = (long) n * n;
        final long rank25 = (count + 3) / 4;
        final long rank75 = (3 * count + 3) / 4;
        final int k25 = pairRank(n, rank25, belowOne, equalOne);
        final int k75 = pairRank(n, rank75, belowOne, equalOne);
        if (k25 >= 0 && k75 >= 0) return selectTwo(lower, k25, k75);
        return new double[] { k25 < 0 ? 1.0 : select(lower, k25),
            k75 < 0 ? 1.0 : select(lower, k75) };
    }

    private static int pairRank(int n, long rank,
            long belowOne, long equalOne) {
        if (rank <= 2 * belowOne) return (int) ((rank - 1) / 2);
        if (rank <= 2 * (belowOne + equalOne) + n) return -1;
        return (int) ((rank - n - 1) / 2);
    }

    /** Partition once while the two requested ranks share a search interval. */
    private static double[] selectTwo(double[] values, int k1, int k2) {
        if (!reuseWorkers) return new double[] { select(values, k1), select(values, k2) };
        int lo = 0, hi = values.length - 1;
        long work = 0, limit = 8L * values.length;
        while (lo < hi) {
            work += hi - lo + 1L;
            if (work > limit) {
                Arrays.sort(values, lo, hi + 1);
                return new double[] {values[k1], values[k2]};
            }
            final double pivot = values[lo + (hi - lo) / 2];
            int lt = lo, cur = lo, gt = hi;
            while (cur <= gt) {
                if (values[cur] < pivot) swap(values, lt++, cur++);
                else if (values[cur] > pivot) swap(values, cur, gt--);
                else cur++;
            }
            if (k2 < lt) hi = lt - 1;
            else if (k1 > gt) lo = gt + 1;
            else {
                final double q1 = k1 < lt ? select(values, k1, lo, lt - 1)
                    : k1 > gt ? select(values, k1, gt + 1, hi) : values[k1];
                final double q2 = k2 < lt ? select(values, k2, lo, lt - 1)
                    : k2 > gt ? select(values, k2, gt + 1, hi) : values[k2];
                return new double[] { q1, q2 };
            }
        }
        return new double[] { values[lo], values[lo] };
    }

    /** In-place three-way quickselect; equal scores need no repeated scans. */
    private static double select(double[] values, int k) {
        return select(values, k, 0, values.length - 1);
    }

    private static double select(double[] values, int k, int lo, int hi) {
        long work = 0, limit = 8L * (hi - lo + 1L);
        while (lo < hi) {
            work += hi - lo + 1L;
            if (work > limit) {
                Arrays.sort(values, lo, hi + 1);
                return values[k];
            }
            final double pivot = values[lo + (hi - lo) / 2];
            int lt = lo, cur = lo, gt = hi;
            while (cur <= gt) {
                if (values[cur] < pivot) swap(values, lt++, cur++);
                else if (values[cur] > pivot) swap(values, cur, gt--);
                else cur++;
            }
            if (k < lt) hi = lt - 1;
            else if (k > gt) lo = gt + 1;
            else return values[k];
        }
        return values[lo];
    }

    private static void swap(double[] values, int i, int j) {
        final double tmp = values[i];
        values[i] = values[j];
        values[j] = tmp;
    }
}
