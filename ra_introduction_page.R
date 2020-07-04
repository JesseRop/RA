tabPanel("About",
         tags$head(tags$script(HTML('
        var fakeClick = function(tabName) {
          var dropdownList = document.getElementsByTagName("a");
          for (var i = 0; i < dropdownList.length; i++) {
            var link = dropdownList[i];
            if(link.getAttribute("data-value") == tabName) {
              link.click();
            };
          }
        };
      '))),
         mainPanel(width = 12,
                   HTML("
<h1>Welcome</h1>
<p>
Here, you can explore our single-cell RNA sequencing data to investigate the phenotypic spectrum of synovial tissue macrophage during development and resolution of rheumatoid arthritis (RA). 
Analyses were performed using the Seurat pipeline and custom scripts which are available upon request. In this web application we use a random subset of 8000 cells for faster computation. Further details can be found in our paper summarized below which you should cite if you use our data.</p>
<p>

</p>

<div class='citation'>
<h4><a href='https://www.google.com' target='_blank'>Distinct synovial tissue macrophage subsets regulate inflammation and remission in
rheumatoid arthritis </a></h4>
<p class='citation-authors'>Stefano Alivernini, Lucy MacDonald, Aziza Elmesmari, Samuel Finlay, Barbara
Tolusso, Maria Rita Gigante, Luca Petricca, Clara Di Mario, Laura Bui, Simone Perniola,
Moustafa Attar, Marco Gessi, Anna Laura Fedele, Sabarinadh Chilaka, Domenico Somma,
Stephen Sansom, Andrew Filer, Charles McSharry, Neal L. Millar, Kristina Kirschner,
Alessandra Nerviani, Myles J. Lewis, Costantino Pitzalis, Andrew R. Clark, Gianfranco
Ferraccioli, Irina Udalova, Christopher D. Buckley, Elisa Gremese, Iain B. McInnes,
Thomas D. Otto and Mariola Kurowska-Stolarska
</p>
</div>


<div>

<h4>Abstract</h4>
<p>Immune-regulatory mechanisms of drug-free remission in rheumatoid arthritis (RA) are unknown.
We hypothesised that synovial tissue macrophages (STM), which persist in remission, contribute
to joint homeostasis. We used single-cell transcriptomics to profile 32000 STMs and identified
phenotypic changes in patients with early/active RA, treatment-refractory/active RA and RA in
sustained remission. Each clinical state was characterised by different frequencies of 9 discrete
phenotypic clusters within 4 distinct STM subpopulations with diverse homeostatic, regulatory and
inflammatory functions. This cellular atlas combined with deep-phenotypic, spatial and functional
analyses of synovial biopsy FACS-sorted STMs revealed two STM subpopulations
(MerTK<sup>pos</sup>/TREM2<sup>high</sup> and MerTK<sup>pos</sup>/LYVE1<sup>pos</sup>) with unique remission transcriptomic signatures
enriched in negative-regulators of inflammation. These STMs were potent producers of
inflammation-resolving lipid mediators and induced the repair response of synovial fibroblasts <i>in
vitro</i>. A low proportion of MerTK<sup>pos</sup> STMs in remission was predictive of flare after treatment
cessation. Therapeutic fostering of MerTK<sup>pos</sup> STM-subpopulations could therefore be a successful
strategy for RA treatment.</p>
</div>

<div class='citation'>
<p> The results are displayed in 2 sections
<ul>
<li><a onclick = 'fakeClick(\"Cluster adjustment\")'> Clustering</a> - Allows cluster visualization and exploration of top cluster markers.</li>
<li><a onclick = 'fakeClick(\"Differential expression\")'> Differential expression (DE)</a> - Comparison of gene expression in macrophages from healthy, UPA, Naive RA, Resistant RA and Remission RA using
parameters specified in our <a href='https://www.google.com' target='_blank'> paper </a> (i.e the first 12 Principal components and 0.5 as the clustering resolution). </li>
</ul>
</p> 
</div>

<div style='border:thin solid black ; padding:30px ; margin:30px'>
<figure>  

<img src='ra_scrna_multiplot_final.png' alt='Results' style='width:100%;'>
<figcaption>scRNAseq defines heterogeneity within MerTK<sup>pos</sup> /CD206<sup>pos</sup> and
MerTK neg /CD206 neg STM populations. <b>(a)</b> UMAP of 9 STM clusters identified by scRNAseq
analysis. <b>(a-h)</b> show data from Healthy (n=4), UPA (n=4), naïve-active RA (n=5), treatment-
resistant RA (n=6) and RA in remission (n=6) in 5 independent experiments. <b>(b)</b> Heatmap of
the top 20 differentially expressed genes per cluster. Top cluster markers and the total
number of genes characterized each cluster are provided. <b>(c)</b> Violin plots represent log-
normalized expression values of STM clusters’ markers with median marked by black dots
while cluster identity by unique colour. <b>(d)</b> Relationship between clusters embedded in the
top 3 Diffusion-map Components. <b>(e)</b> Hierarchical clustering of STMs. <b>(f)</b> MerTK expression
in the 9 STM clusters. <b>(g)</b> Proposed classification of human STMs. <b>(h)</b> Split UMAP and dot
plots of relative changes in the STM clusters between groups. Significant differences
(*p<0.05) between the given condition and at least two other conditions in Two-way ANOVA
with Tukey’s correction. Precise p-values in Supplementary Figure 3a which can be found in our <a href='https://www.google.com' target='_blank'> paper </a>.</figcaption>
</figure>
</div>

<footer align = 'center' style = 'position:absolute; bottom:80; width:95%; height:20px; padding: 10px; z-index: 200; background-color:PowderBlue;'> <hr></hr> <p>This web-app (version 1.01) was developed by <a href='https://https://www.gla.ac.uk/researchinstitutes/iii/staff/?webapp=staffcontact&action=person&id=4edfe8ec8a95' target='_blank'> Jesse Rop </a> under the supervision of <a href='https://https://www.gla.ac.uk/researchinstitutes/iii/staff/thomasdanotto/' target='_blank'> Dr Thomas Otto </a> and with the assistance of <a href='https://www.linkedin.com/in/lucy-macdonald-07248b177/?originalSubdomain=uk' target='_blank'> Lucy Macdonald </a> and <a href='https://https://www.gla.ac.uk/researchinstitutes/iii/staff/?webapp=staffcontact&action=person&id=49d8ebe38396' target='_blank'> Scott Arkinson </a>. Further acknowledgements on the generation of the data hosted in the app can be found in our <a href='https://www.google.com' target='_blank'> paper </a>.</p> <p> The code used to develop the app can be found <a href='https://github.com/JesseRop/' target='_blank'> here </a>. For any queries or suggestions kindly <a href='mailto:ThomasDan.Otto@glasgow.ac.uk' target='_blank'> contact us</a>. An updated version coming soon.</p>
</footer>"))
         
)
