# Computer Vision (AIT-SP17)
Author: Sooji Son

## File Structure

```
.
├── README.md
├── .gitignore
├── runme.bat - a Windows script that you will fill in to demonstrate execution of your program
├── runme.sh - same as <code>runme.bat, but for Mac OS X
└── src/ - direcotry with source code
    ├── Makefile
    ├── imagepro.[vcproj/sln/suo] - for Visual Studio 2005 on Windows
    ├── imgpro.cpp - Main, parses the command line args, calls the appropriate image function
    ├── R2Image.[cpp/h] - Image class with processing functions (to be edited by me)
    ├── R2Pixel.[cpp/h] - Pixel class
    ├── R2/ - A library of useful 2D geometric primitives
    └── jpeg/ - A library for reading/writing JPEG files
```

## Compilation
If you are developing on a Windows machine and have Visual Studio
installed, use the provided project solution file (assn1.sln) in the
src/ directory to build the program. If you are developing on a Mac or
Linux machine, cd into the src/ directory and type "make". In either
case, an executable called imgpro (or imgpro.exe) will be created in
the src/ directory.
