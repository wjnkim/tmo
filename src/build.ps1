param(
    [string]$StataDirectory = "",
    [string]$JavaHome = "",
    [string]$ScratchDirectory = ""
)
$ErrorActionPreference = "Stop"
if (!$StataDirectory) { $StataDirectory = $env:STATA_HOME }
if (!$StataDirectory) {
    $cpdStataInstalls = @(Get-ChildItem -Path (Join-Path $env:ProgramFiles 'Stata*') -Directory -ErrorAction SilentlyContinue |
        Where-Object { Test-Path -LiteralPath (Join-Path $_.FullName 'utilities/jar/sfi-api.jar') } | Sort-Object FullName -Descending)
    if ($cpdStataInstalls.Count) { $StataDirectory = $cpdStataInstalls[0].FullName }
}
if (!$StataDirectory) { throw 'Stata was not found. Set STATA_HOME or supply -StataDirectory.' }
if (!$JavaHome) { $JavaHome = $env:JAVA_HOME }
if (!$JavaHome) {
    $cpdJdks = @(Get-ChildItem -Path (Join-Path $StataDirectory 'utilities/java/windows-x64/*') -Directory -ErrorAction SilentlyContinue |
        Where-Object { Test-Path -LiteralPath (Join-Path $_.FullName 'bin/javac.exe') } | Sort-Object FullName -Descending)
    if ($cpdJdks.Count) { $JavaHome = $cpdJdks[0].FullName }
}
if (!$JavaHome) { throw 'JDK 11 or newer is required. Set JAVA_HOME or supply -JavaHome.' }
$cpdCompiler = Join-Path $JavaHome "bin\javac.exe"
$cpdJarTool = Join-Path $JavaHome "bin\jar.exe"
$cpdApi = Join-Path $StataDirectory "utilities\jar\sfi-api.jar"
$cpdScratch = if ($ScratchDirectory) {
    [IO.Path]::GetFullPath($ScratchDirectory)
} else {
    [IO.Path]::GetFullPath([IO.Path]::GetTempPath())
}
if (!(Test-Path -LiteralPath $cpdScratch -PathType Container)) {
    throw "Provide an existing -ScratchDirectory for temporary build files."
}
$cpdBuild = Join-Path $cpdScratch ("tmo_java_build_" + [Guid]::NewGuid().ToString("N"))
$cpdClasses = Join-Path $cpdBuild "classes"
$cpdBuiltJar = Join-Path $cpdBuild "tmo.jar"
try {
    New-Item -ItemType Directory -Path $cpdClasses -Force | Out-Null
    & $cpdCompiler --release 11 -cp $cpdApi -d $cpdClasses (Join-Path $PSScriptRoot "tmo.java")
    if ($LASTEXITCODE -ne 0) { throw "tmo.java compilation failed." }
    & $cpdJarTool --create --file $cpdBuiltJar -C $cpdClasses tmo
    if ($LASTEXITCODE -ne 0) { throw "tmo.jar creation failed." }
    Copy-Item -LiteralPath $cpdBuiltJar -Destination (Join-Path $PSScriptRoot "tmo.jar") -Force
    Write-Output "Built src/tmo.jar. Use a fresh Stata process to load the rebuilt class."
} finally {
    $cpdResolvedBuild = [IO.Path]::GetFullPath($cpdBuild)
    $cpdScratchPrefix = $cpdScratch.TrimEnd('\') + '\'
    if (!$cpdResolvedBuild.StartsWith($cpdScratchPrefix, [StringComparison]::OrdinalIgnoreCase) -or
            [IO.Path]::GetFileName($cpdResolvedBuild) -notlike "tmo_java_build_*") {
        throw "Refusing cleanup outside the verified temporary build directory."
    }
    if (Test-Path -LiteralPath $cpdResolvedBuild) {
        Remove-Item -LiteralPath $cpdResolvedBuild -Recurse -Force
    }
}
