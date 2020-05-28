
getName <- function(v) {
deparse(substitute(v))
}

## Draw a small histogram of each different dataset, separated by animal/field of view

dir.create("animal")
dir.create("fov")

# Plot histograms of log2 H2BGFP intensities (all data; separated per animal)
for (dataset in names(everything)) {
		name <- sprintf("animal/%s.png",dataset)
		print(sprintf("Generating %s", name))
		title <- sprintf("Histone dilution in dataset %s",dataset)
		png(name)
		hist(log2(everything[[dataset]]),main=title,xlab="log2(Histone Intensity)")
		dev.off()
}

# Plot histograms of log2 H2BGFP intensities (all data; separated per field-of-view)
for (dataset in names(fov)) {
		name <- sprintf("fov/%s.png",dataset)
		print(sprintf("Generating %s", name))
		title <- sprintf("Histone dilution in dataset %s",dataset)
		png(name)
		hist(log2(fov[[dataset]]),main=title,xlab="log2(Histone Intensity)")
		dev.off()
}
