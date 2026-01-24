# Vector Class Library

## Required: Agner Fog's Vector Class Library

RECURSUM requires the Vector Class Library for SIMD support.

### Installation Options

**Option 1: Download Manually**

Download from: https://github.com/vectorclass/version2

Copy `vectorclass.h` and related files to this directory:
```
include/recursum/vectorclass.h
include/recursum/vectorclass/...
```

**Option 2: Git Submodule**

```bash
cd RECURSUM
git submodule add https://github.com/vectorclass/version2 external/vectorclass
ln -s ../../external/vectorclass/vectorclass.h include/recursum/vectorclass.h
```

**Option 3: System Installation**

```bash
# If vectorclass is installed system-wide
# CMake will find it automatically via CMAKE_INCLUDE_PATH
```

## License

Vector Class Library is licensed under Apache License 2.0.

RECURSUM uses this library under the terms of the Apache 2.0 license.
See: https://github.com/vectorclass/version2/blob/master/LICENSE

## Attribution

Vector Class Library
Copyright (c) 2012-2021 Agner Fog
https://github.com/vectorclass/version2
