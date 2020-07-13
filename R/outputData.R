

outputData = function(clusters, stories, plotDirectory, genome='hg19', resourceDirectory='superFreqResources') {
    dataDir = paste0(plotDirectory, '/data')
    ensureDirectoryExists(dataDir)
    
    #copy numbers segments. bed with extra column for gene list.
    load(paste0(resourceDirectory, '/COSMIC/cosmicCounts.Rdata'))
	cosmicGenes = names(cosmicCounts$geneCounts)
    lapply(names(clusters), function(sample) {
        file = paste0(dataDir, '/CNAsegments_', sample, '.tsv')
        #if ( file.exists(file) ) return()
        cna = clusters[[sample]]
        genesBySegment = lapply(1:nrow(cna$clusters), function(i) rownames(cna$CR)[cna$CR$x1 <= cna$clusters$x2[i] & cna$CR$x2 >= cna$clusters$x1[i]])
        cosmicGenesBySegment = lapply(genesBySegment, function(genes) genes[genes %in% cosmicGenes])
        geneText = sapply(genesBySegment, function(genes) gsub(' ', ',', do.call(paste, as.list(genes))))
        cosmicGeneText = sapply(cosmicGenesBySegment, function(genes) {
        	if ( length(genes) == 0 ) return('')
        	return(gsub(' ', ',', do.call(paste, as.list(genes))))
        })
        chr = xToChr(cna$clusters$x1, genome=genome)
        start = xToPos(cna$clusters$x1, genome=genome)
        end = xToPos(cna$clusters$x2, genome=genome)
        df = cbind(chr, start, end, cna$clusters, genes=geneText, COSMIC_genes=cosmicGeneText)

        write.table(df, file=file, sep='\t', quote=F, row.names=F)
    })
        
    
    #LFC and BAF by gene (sorted by position)
    lapply(names(clusters), function(sample) {
        file = paste0(dataDir, '/CNAbyGene_', sample, '.tsv')
        if ( file.exists(file) ) return()
        cna = clusters[[sample]]
        chr = xToChr(cna$CR$x1, genome=genome)
        start = xToPos(cna$CR$x1, genome=genome)
        end = xToPos(cna$CR$x2, genome=genome)
        gene = rownames(cna$CR)
        df = cbind(chr, start, end, gene, cna$CR)

        write.table(df, file=file, sep='\t', quote=F, row.names=F)
    })

    #stories from clonal tracking, ID, list of immediate subclones, clonality, error
    lapply(names(stories), function(ind) {
        file = paste0(dataDir, '/clones_', ind, '.tsv')
        if ( file.exists(file) ) return()
        story = stories[[ind]]
        df = cbind('ID'=rownames(story$consistentClusters$cloneStories), story$consistentClusters$cloneStories[,c('stories', 'errors')])

        write.table(df, file=file, sep='\t', quote=F, row.names=F)
    })
    
    

}
