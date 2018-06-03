## Changelog

### v2.1.0

- zopflipng: scaled down input file size ranges
- zopflipng: scaled down number of iterations in file size ranges

### v2.0.0

- added more granular definitions of number of zopfli iterations based upon file size (defined in 9 size ranges).  These ranges were defined by experimentation with quantized PNG files that were created with pngquant
- removed `-m` option
- added non-zero exit status codes for zopflipng command line errors
- added gitignore file
- added zopflipng version reporting in help text and with `-v` + `--version` flags
- minor help text updates

### Fork

This repository was forked from the google/zopfli repository at commit https://github.com/chrissimpkins/zopfli/commit/ae43a8b73827577c4b19b005b6eed81f5cf9bbac (version tag zopfli-1.0.2)