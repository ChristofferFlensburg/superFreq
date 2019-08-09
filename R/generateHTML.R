
#' generates an HTML to navigate superFreq output
#'
#' @param metaDataFile Character. A path to the meta data file used for the run. 
#' @param outputFile Character. A path to where the HTML is stored. Has to be inside the plotDirectory for now.
#'
#' @details This function is in pre-alpha development as of version 0.9.23. Feel free (please!) to add to this on github even without any html experience. I dont have any. The function generates a .html that should link the plots and explanation of how to interpret the plots, how they were generated and how to get the raw data for further analysis. Ideally there should be a good landing page as well for overview.
#'          
#'
#' @export
#'
#' @examples
#' \dontrun{
#' metaDataFile = 'metaData.tsv'
#' plotDirectory='plots'
#' 
#' printHTML(metaDataFile=metaDataFile, outputFile=paste0(plotDirectory, '/superFreq.html'))
#' }
printHTML = function(metaDataFile, outputFile) {
  metaData = importSampleMetaData(metaDataFile)
  script = generateHtml(metaData)
  cat(script, file=outputFile)
}

generateHtml = function(metaData) {
  head = generateHead()
  body = generateBody(metaData)


  script = pasteLines('<!DOCTYPE html>',
    '<html>', '',
    head, '',
    body,
    '</body>',
    '</html>\n')

  return(script)
}









generateTabs =function() {
  tabs = '<div class="tab">
  <button class="tablinks" onclick="openMain(event, \'cna\')" id="defaultOpen">CNA</button>
  <button class="tablinks" onclick="openMain(event, \'scatters\')">SNV scatters</button>
  <button class="tablinks" onclick="openMain(event, \'river\')">River</button>
</div>
'
  return(tabs)
}




generateBody = function(metaData) {
  tabs = generateTabs()

  overviewTab = generateOverviewTab(metaData)
  cnaTab = generateCnaTab(metaData)
  scatterTab = generateScatterTab(metaData)
  riverTab = generateRiverTab(metaData)

  script = generateScript()

  body = pasteLines(tabs, '',
    cnaTab, '',
    scatterTab, '',
    riverTab, '',
    script, '')

  return(body)
}

generateOverviewTab = function(metaData) {
  preambl = '<div id="cna" class="tabcontent">
<h3>Clonal structure and key mutations</h3>

<h3>Top somatic SNVs (fold out table, 5-10 hits)</h3>

<h3>Copy number profiles by sample (fold out, only top panel)</h3>
'

  overviewTab = pasteLines(preambl, '</div>\n')

  return(overviewTab)
}


generateCnaTab = function(metaData) {
  preambl = '<div id="cna" class="tabcontent">
<div data-role="collapsible">
<h1>Help</h1>
<div data-role="collapsible">
<h2>Reading the CNA plots</h2>
<p>These plots show the copy number calls of each sample.</p>
<p>The top panel shows the coverage log fold change (LFC) with respect to the reference normal samples. The small dots represent the LFC of a gene, and size of the dots represents the condifence of the LFC, as estimated by limma-voom. The segmentation is shown by the opaque dots and lines, with the horizontal lines showing the range of the segment, and the size of dot and error bar showing the confidence in the LFC. This is normalised for ploidy so that each segment fit as well as possible with some copy number call of some clonality. The SNP frequencies are also used in the ploidy estimate.</p>
</div>

<div data-role="collapsible">
<h2>How superFreq calls CNAs</h2>
<p>Like most modern copy number callers, SuperFreq uses both read depth and allele frequency of germline heteroygous SNPs to call copy numbers. The read depth is analysed by gene with limma-voom, a state of the art differential expression method normally used for RNA-seq, to get good estimates of read depth log fold change compared to the reference normals. So each sample is compared one-to-many against the reference normals. The heterozygous germline SNPs are identified from the matched normal if available, or from the sample itself otherwise (restrciting to VAF between 5% and 95%). LINK TO READ DEPTH QC AND VAF QC PLOTS.</p>

<p>The genes then undergo hierarchical clustering by chromosome, where the most similar (within error estimates) neighbours, in both LFC and MAF are paired up until there are no more sufficiently similar neighbouring segments.</p>

<p>Superfreq allows each segment to separate clonalities. THe ploidy is decided so that each segments fit a copy number call at some clonality for both LFC and MAF, within uncertainty. This process is shown in the LINK TO MAYPOLE PLOT. The allele sensitive copy number and clonality is then decided for each segment from the best fit in the maypole plot.</p>
</div>

<div data-role="collapsible">
<h2>Getting the data</h2>
<p>The log fold changes are available in <a href="differentRegionsSamples.xls">
this file</a>, or by exon in <a href="differentRegionsSamplesExon.xls"> this file</a>. The output from limma-voom is in RDATA FIT. The normalised LFC and MAF by gene are available in RDATA CLUSTER.</p>
</div>
</div>
'
  imgs = lapply(metaData$NAME, function(sample) paste0('<h4>', sample, '</h4>\n<img src="CNV/', sample, '.png" alt="a CNV plot" title="coverage log fold change, minor allele frequency of heterozygous germline SNPs, clonality of CNA call" width="800px" height="400px">\n<br>'))
  #add collapsible panel with the indivdual chromosomes.
  imgs = do.call(pasteLines, imgs)
  cnaTab = pasteLines(preambl, imgs, '</div>\n')

  return(cnaTab)
}

generateScatterTab = function(metaData) {
  preambl = '<div id="scatters" class="tabcontent">
<div data-role="collapsible">
<h1>Help</h1>
<div data-role="collapsible">
<h2>Reading the scatter plots</h2>
<p>talk about classes of variants and how they show up in the scatters. Mayb mention contamination?</p>
</div>

<div data-role="collapsible">
<h2>How superFreq filters variants and calls somatic variants</h2>
<p>asdf</p>
</div>

<div data-role="collapsible">
<h2>Getting the data</h2>
<p>point to somaticVariants, and .Rdata.</p>
</div>
</div>
'
  samples = metaData$NAME
  imgs = lapply(1:(length(samples)-1), function(i) {
    lapply((i+1):length(samples), function(j) paste0('<img src="scatters/', samples[i], '/', samples[j], '/all.png" alt="scatter plot" title="scatter plot of variant allele frequencies of SNVs and short indels" width="500px" height="500px">\n<br>'))
  })
  imgs = do.call(pasteLines, as.list(unlist(imgs)))
  scatterTab = pasteLines(preambl, imgs, '</div>\n')

  return(scatterTab)
}

generateRiverTab = function(metaData) {
  preambl = '<div id="river" class="tabcontent">\n<h3>River</h3>
<div data-role="collapsible">
<h1>Help</h1>
<div data-role="collapsible">
<h2>Reading the river plots</h2>
<p>Talk about visualisations of clones. Maybe an example with arrows etc? Link to some powerpoint presentation? Talka bout filtering and dodgy rivers.</p>
</div>

<div data-role="collapsible">
<h2>How superFreq tracks clones</h2>
<p>difference between dodgy and non-dodgy.</p>
</div>

<div data-role="collapsible">
<h2>Getting the data</h2>
<p>point to .xls in rivers and .Rdata.</p>
</div>
</div>
'
  imgs = lapply(metaData$NAME, function(sample) paste0('<img src="CNV/', sample, '.png" alt="a CNV plot" title="coverage log fold change, minor allele frequency of heterozygous germline SNPs, clonality of CNA call">'))
  imgs = do.call(pasteLines, imgs)
  cnaTab = pasteLines(preambl, imgs, '</div>\n')

  return(cnaTab)
}


pasteLines = function(...) paste(..., sep='\n')


generateScript = function() {
  script = '<script>
function openMain(evt, tabName) {
    var i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = "none";
    }
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(" active", "");
    }
    document.getElementById(tabName).style.display = "block";
    evt.currentTarget.className += " active";
}

// Get the element with id="defaultOpen" and click on it
document.getElementById("defaultOpen").click();
</script>
'

  return(script)
}

generateHead = function() {
  head = '<head>
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" href="https://code.jquery.com/mobile/1.4.5/jquery.mobile-1.4.5.min.css">
<script src="https://code.jquery.com/jquery-1.11.3.min.js"></script>
<script src="https://code.jquery.com/mobile/1.4.5/jquery.mobile-1.4.5.min.js"></script>

<style>
* {box-sizing: border-box}
body {font-family: "Lato", sans-serif;}

/* Style the tab */
div.tab {
    float: left;
    border: 1px solid #ccc;
    background-color: #f1f1f1;
    width: 30%;
    height: 300px;
}

/* Style the buttons inside the tab */
div.tab button {
    display: block;
    background-color: inherit;
    color: black;
    padding: 22px 16px;
    width: 100%;
    border: none;
    outline: none;
    text-align: left;
    cursor: pointer;
    transition: 0.3s;
    font-size: 17px;
}

/* Change background color of buttons on hover */
div.tab button:hover {
    background-color: #ddd;
}

/* Create an active/current "tab button" class */
div.tab button.active {
    background-color: #ccc;
}

/* Style the tab content */
.tabcontent {
    float: left;
    padding: 0px 12px;
    border: 1px solid #ccc;
    width: 70%;
    border-left: none;
}
</style>
</head>
'
  return(head)
}
