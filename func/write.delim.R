# Function:
#	write.delim: writes tab-delimited text file
# Arguments:
#	x, file, row.names, col.names, quote: as in utils::write.table

write.delim <- function(x, file, row.names=FALSE, col.names=TRUE, quote=FALSE) {
	write.table(x, file, row.names=row.names, col.names=col.names, quote=quote, sep="\t")
}
