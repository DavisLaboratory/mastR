# mastR 0.99.0

* Added a `NEWS.md` file to track changes to the package.

# mastR 0.99.1

* Fixed notes in R CMD check.

# mastR 0.99.2

* Added non-emtpy return value to man files and update pca plot function.

# mastR 0.99.3

* Added *.Rproj into .gitignore.

# mastR 0.99.4

* Deleted *.Rproj file to stop git track.

# mastR 0.99.5

* Updated R dependency, delete .Rhistory file, move all generic functions to a separate file "AllGenerics.R".

# mastR 0.99.6

* Removed paste in message().

# mastR 0.99.7

* Improved test coverage depth, modified demo vignette and update Suggests packages.

# mastR 0.99.8

* Linked external packages in vignette using `BiocStyle` functions, fixed the note of using `paste` in conditions.

# mastR 0.99.9

* Updated vignette for `BiocStyle` packages link functions.

# mastR 1.2.1

* Specify Matrix (<= 1.6.1.1) to avoid conflicts between SeuratObject and Matrix 1.6-2. Will fix to update to the latest version in BiocVersion 3.19.

# mastR 1.3.1

* Remove Matrix version bound.

# mastR 1.3.2

* Update function `voom_fit_treat()` to allow pass user-defined contrast matrix.

# mastR 1.3.3

* Fix names mismatch problem when passing user-defined contrast matrix to DE analysis functions.

# mastR 1.3.4

* Update function `sig_gseaplot()` to allow more custom arguments for `enrichplot::gseaplot2()`.

# mastR 1.3.5

* Update function `remove_bg_exp()` using Gaussian distribution percentiles to replace min-max scaling as relative exppression within each sample.

# mastR 1.3.6

* Add bioRxiv citation for mastR package.

# mastR 1.9.1

* Add Bioinformatics citation for mastR package, add sticker and remove BisqueRNA dependency.
