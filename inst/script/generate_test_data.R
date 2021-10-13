## making test data
set.seed(10)

test <- readRDS("mouse_m6A.rds")
GenomicRanges::elementMetadata(test) <- NULL
test <- test[sample(seq_len(length(test)),100)]
test <- test+10
export.gff3(test, "inst/extdata/test.rds")

## generate the test gff3
anno <- rtracklayer::import.gff3("gencode.vM23.primary_assembly.annotation.gff3")
anno_exon <- anno[anno$type == "exon"]
anno_exon <- as.data.frame(anno_exon) %>% group_by(transcript_id) %>%
  mutate(trans_len = sum(width)) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
anno$trans_len <- anno_exon$trans_len[
  match(anno$transcript_id, anno_exon$transcript_id)]

anno_exon <- anno_exon[order(anno_exon$level,
                             anno_exon$transcript_support_level,
                             -anno_exon$trans_len)]
k <- findOverlaps(test, anno_exon, select = "first")
k2 <- anno_exon[k]

k3 <- anno[anno$transcript_id %in% c(k2$transcript_id,k22$transcript_id)]

export.gff3(k3, "inst/extdata/annotation_test.gff3")