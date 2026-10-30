# TMO Java kernel

Normal Stata installation uses the included tmo.jar. The public tmo.ado command falls back to its Mata implementation when the Java kernel is unavailable. Users do not need a JDK or a source build.

## Build from source on Windows

The PowerShell build script compiles tmo.java against Stata's sfi-api.jar, targets Java 11 bytecode, and replaces src/tmo.jar. From the repository root, run:

    ./src/build.ps1

The script detects a standard Windows Stata installation and its bundled JDK. Set STATA_HOME or JAVA_HOME, or pass -StataDirectory and -JavaHome, for other installations. The script uses the system temporary directory by default; -ScratchDirectory accepts an existing alternative directory. Temporary classes are removed after the build. Start a fresh Stata process before loading a rebuilt jar.

## Check the kernel

The source checkout includes deterministic source and shipped-jar checks:

    ./src/checks/run.ps1 -Mode Source
    ./src/checks/run.ps1 -Mode Jar

Both modes use the system temporary directory by default. Pass -ScratchDirectory to use another existing directory and -JavaHome for a different JDK. The checks compile only temporary test classes and do not replace tmo.jar. Public-command Stata checks are also needed to verify integration with Stata estimators and the Stata Function Interface.
