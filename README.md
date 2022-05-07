## zopfli

**This repository is a derivative of the [google/zopfli source repository](https://github.com/google/zopfli) and includes modifications of the upstream source files in the master branch.**

Zopfli Compression Algorithm is a compression library programmed in C to perform
very good, but slow, deflate or zlib compression.

The basic function to compress data is ZopfliCompress in zopfli.h. Use the
ZopfliOptions object to set parameters that affect the speed and compression.
Use the ZopfliInitOptions function to place the default values in the
ZopfliOptions first.

ZopfliCompress supports deflate, gzip and zlib output format with a parameter.
To support only one individual format, you can instead use ZopfliDeflate in
deflate.h, ZopfliZlibCompress in zlib_container.h or ZopfliGzipCompress in
gzip_container.h.

 - ZopfliDeflate creates a valid deflate stream in memory, see:
http://www.ietf.org/rfc/rfc1951.txt
 - ZopfliZlibCompress creates a valid zlib stream in memory, see:
http://www.ietf.org/rfc/rfc1950.txt
 - ZopfliGzipCompress creates a valid gzip stream in memory, see:
http://www.ietf.org/rfc/rfc1952.txt

This library can only compress, not decompress. Existing zlib or deflate
libraries can decompress the data.

zopfli_bin.c is separate from the library and contains an example program to
create very well compressed gzip files. Currently the makefile builds this
program with the library statically linked in.

The source code of Zopfli is under src/zopfli. Build instructions:

To build zopfli, compile all .c source files under src/zopfli to a single binary
with C, and link to the standard C math library, e.g.:
`gcc src/zopfli/*.c -O2 -W -Wall -Wextra -Wno-unused-function -ansi -pedantic -lm -o zopfli`

A makefile is provided as well, but only for linux. Use `make` to build the
binary, `make libzopfli` to build it as a shared library. For other platforms,
please use the build instructions above instead.

Zopfli Compression Algorithm was created by Lode Vandevenne and Jyrki
Alakuijala, based on an algorithm by Jyrki Alakuijala.

## zopflipng

ZopfliPNG is a command line program to optimize the Portable Network Graphics
(PNG) images. This version has the following features:

- uses Zopfli compression for Deflate compression,
- compares several strategies for choosing scanline filter codes,
- chooses a suitable color type to losslessly encode the image,
- removes all chunks that are unimportant for the typical web use (metadata,
  text, etc...),
- optionally alters the hidden colors of fully transparent pixels for more
  compression, and,
- optionally converts 16-bit color channels to 8-bit.

Important:

This PNG optimizer removes ancillary chunks (pieces of metadata) from the
PNG image that normally do not affect rendering. However in special
circumstances you may wish to keep some. For example for a design using
custom gamma correction, keeping it may be desired. Visually check in the
target renderer after using ZopfliPNG. Use --keepchunks to keep chunks, e.g.
--keepchunks=gAMA,pHYs to keep gamma and DPI information. This will increase
file size. The following page contains a list of ancillary PNG chunks:
http://www.libpng.org/pub/png/spec/1.2/PNG-Chunks.html

Build instructions:

To build ZopfliPNG, compile all .c, .cc and .cpp files from src/zopfli,
src/zopflipng and src/zopflipng/lodepng, except src/zopfli/zopfli_bin.c, to a
single binary with C++, e.g.:
`g++ src/zopfli/{blocksplitter,cache,deflate,gzip_container,hash,katajainen,lz77,squeeze,tree,util,zlib_container,zopfli_lib}.c src/zopflipng/*.cc src/zopflipng/lodepng/*.cpp -O2 -W -Wall -Wextra -Wno-unused-function -ansi -pedantic -o zopflipng`

A makefile is provided as well, but only for linux: use `make zopflipng` with
the Zopfli makefile. For other platforms, please use the build instructions
above instead.

The main compression algorithm in ZopfliPNG is ported from WebP lossless, but
naturally cannot give as much compression gain for PNGs as it does for a more
modern compression codec like WebP
https://developers.google.com/speed/webp/docs/webp_lossless_bitstream_specification.

Compared to libpng -- an often used PNG encoder implementation -- ZopfliPNG uses
2-3 orders of magnitude more CPU time for compression. Initial testing using a
corpus of 1000 PNGs with translucency, randomly selected from the internet,
gives a compression improvement of 12% compared to convert -q 95, but only 0.5%
compared to pngout (from better of /f0 and /f5 runs).

By releasing this software we hope to make images on the web load faster without
a new image format, but the opportunities for optimization within PNG are
limited. When targeting Android, Chrome, Opera, and Yandex browsers, or by using
suitable plugins for other browsers, it is good to note that WebP lossless
images are still 26 % smaller than images recompressed with ZopfliPNG.
