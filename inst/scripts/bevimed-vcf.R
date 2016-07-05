get_sample_names <- function(file, connection_type=base::file) {
	cnx <- connection_type(file, open="r")
	on.exit(close(cnx))
	headings <- NULL
	while (is.null(headings)) {
		l <- readLines(cnx, n=1)
		if (grepl(x=l, pattern="^#CHROM")) headings <- strsplit(l, split="\t")[[1]]
	}
	headings[10:length(headings)]
}

compressed_vcf_sample_names <- function(vcf_file_name) {
	get_sample_names(paste0("zcat ", vcf_file_name), connection_type=pipe)
}

just_counts <- function(parts, description_columns=9) {
	y <- sapply(parts, "[", -(1:description_columns))
	structure(grepl(x=y, pattern="^[^0.][/|].") + grepl(x=y, pattern="^.[/|][^0.]"), dim=c(if (length(parts) > 0) length(parts[[1]])-description_columns else 0, length(parts)))
}

var_info <- function(parts, description_columns=9) structure(dimnames=list(NULL, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","FORMAT","INFO")), structure(dim=c(length(parts),description_columns), t(sapply(parts, "[", 1:description_columns))))

get_block_parts <- function(vcf_file_name, chr, from, to) {
	cmd <- paste("tabix ", vcf_file_name, " ", chr, ":", from, "-", to, sep="")
	z <- pipe(cmd)
	lines <- grep(value=TRUE, pattern="^#", invert=TRUE, x=readLines(z))
	close(z)
	strsplit(lines, split="\t")
}

bevimed_vcf <- function(vcf_file_name, chr, from, to, case_IDs, samples=compressed_vcf_sample_names(vcf_file_name), description_columns=9, min_ac=1, assoc_prior=0.01, filter_variants=function(info) rep(TRUE, nrow(info)), ...) {
	file_samples <- compressed_vcf_sample_names(vcf_file_name) 
	sample_inds <- match(samples, file_samples)
	parts <- get_block_parts(vcf_file_name, chr, from, to)
	counts <- t(just_counts(parts, description_columns))
	info <- var_info(parts, description_columns)
	selected_variants <- filter_variants(info)
	stopifnot(sum(selected_variants) > 0)

	x <- summary(bevimed(
		y=samples %in% case_IDs, 
		G=counts[selected_variants,sample_inds,drop=FALSE], 
		min_ac=min_ac,
		...
	), gamma1_prior=assoc_prior)
	x[["conditional_prob_pathogenic"]] <- setNames(x[["conditional_prob_pathogenic"]], paste0(info[selected_variants,1],":",info[selected_variants,2]))
	x
}

windows_x_phenotypes <- function(vcf_file_name, chr, from, to, list_of_case_sets, min_ac_per_case_set=rep(1, length(list_of_case_sets)), samples=compressed_vcf_sample_names(vcf_file_name), window_names=paste0(chr, ":", from, "-", to), description_columns=9, mc.cores=1L, ...) {
	file_samples <- compressed_vcf_sample_names(vcf_file_name)
	sample_inds <- match(samples, file_samples)
	result_mat <- t(mcmapply(
		FUN=function(chr, from, to) {
			parts <- get_block_parts(vcf_file_name, chr, from, to)
			counts <- t(just_counts(parts, description_columns))
			if (prod(dim(counts)) == 0) rep(NA, length(list_of_case_sets)) 
			else mapply(
				function(cases, min_ac) log_BF(
					y=samples %in% cases, 
					G=counts[,sample_inds,drop=FALSE], 
					min_ac=min_ac,
					...
				),
				list_of_case_sets,
				min_ac_per_case_set
			)
		},
		chr,
		from,
		to,
		mc.cores=mc.cores,
		USE.NAMES=FALSE
	))
	rownames(result_mat) <- window_names
	result_mat
}
