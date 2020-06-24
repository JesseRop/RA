
#setwd("D:/GCRF_UoG/Vicky_JCR_Shiny/Tcells_amp_adopted")
source("ra_libs_raw.R")
##Creating the ui object for the user interface

ui = navbarPage(title = "Synovial atlas of Macrophages in RA", theme = "amp_custom.css", 
                #tags$head(tags$style("div.dataTables_scrollHead span {color: black}")),
                
                source(file.path("ra_introduction_page.R"), local = TRUE)$value,
                
                ##The first tab which is for the differential expression visualization
                tabPanel("Clustering",
                         value = "Cluster adjustment",
                         wellPanel(
                             fluidRow(
                                 
                                 column(
                                     12,
                                     h4("Cluster adjustment"),
                                     tabsetPanel(
                                         #title = "UMAP Clustering",
                                         # The id lets us use input$tabset1 on the server to find the current tab
                                         #id = "vis_de", 
                                         #height = "250px",
                                         #width = 12,
                                         tabPanel("Labelling populations", 
                                                  wellPanel(
                                                      #sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = res1, min = res1, max = res2, step = diff_res, round = F),
                                                      withSpinner(plotOutput("labelled_umap"))
                                                      
                                                      # )
                                                      ,
                                                      tags$hr(),
                                                      uiOutput("cluster_annot")
                                                      # )
                                                      
                                                  )),
                                         tabPanel("Compare groups", withSpinner(plotOutput("all_groups")))
                                         
                                     ),
                                     tags$hr(),
                                     wellPanel(style = "background:#385A4F",
                                               tags$hr(),
                                               tags$p(style = "font-family:Arial;color:white",
                                                      paste("Adjusting resolution (0.15, 0.25, 0.35, 0.45 or 0.55) to set the clustering granularity. This option allows the sub-division of clusters further into sub-populations and the subsequent interrogation of differential expression. Labelling of cell populations can then be done based on top cluster markers in the 'Cluster marker' table on the left.")
                                                      
                                               ),
                                               tags$hr()
                                     )
                                 )
                             ),
                             fluidRow(
                                 column(
                                     12, 
                                     h4("Cluster markers"),
                                     tabsetPanel(
                                         #title = "Cluster markers",
                                         # The id lets us use input$tabset1 on the server to find the current tab
                                         # id = "vis_de", 
                                         # #height = "250px",
                                         # width = 12,
                                         tabPanel("Marker table", 
                                                  wellPanel(
                                                      # box(
                                                      # width = NULL,
                                                      # solidHeader = TRUE,
                                                      uiOutput("dyn_clusters"),
                                                      uiOutput("topclgenes"),
                                                      withSpinner(DTOutput("top_conserved_genes"))
                                                  )
                                                  
                                                  # )
                                         ),
                                         tabPanel("Marker feature plots", 
                                                  uiOutput("top_markers_umap"),
                                                  withSpinner(plotOutput("conserved_markers_umap")))
                                     ),
                                     tags$hr(),
                                     uiOutput("box_2_2"),
                                     HTML(
                                         "<div class='text-center'><button class='btn-success btn-lg' onclick = 'fakeClick(\"Differential expression\")'> Explore differential expression in generated clusters</button></div>"
                                     )
                                 )
                             )
                             
                         )
                ),
                tabPanel("Differential expression (DE)",
                         value = "Differential expression",
                         tabsetPanel(
                             tabPanel("Gene view",
                                      value = "gene_view",
                                      wellPanel(
                                          
                                          fluidRow(
                                              #h3("Gene expression"),
                                              wellPanel(
                                                  fluidRow(
                                                      column(
                                                          12,
                                                          h4("Single gene DE visualization"),
                                                          selectInput(inputId = "de_genes", label = strong("Choose gene:"),choices = all_genes_ra, multiple = F, selected = choice_gene),
                                                      )
                                                  ),
                                                  fluidRow(
                                                      tabsetPanel(
                                                          #title = "Single gene differential expression",
                                                          # The id lets us use input$tabset1 on the server to find the current tab
                                                          #id = "vis_de", 
                                                          #height = "250px",
                                                          #width = 12,
                                                          tabPanel("Violin plot", withSpinner(plotOutput("de_stim_vs_ctrl_vp"))),
                                                          tabPanel("UMAP feature plot", withSpinner(plotOutput("de_stim_vs_ctrl_um")))
                                                      ),
                                                      
                                                      
                                                  ),
                                                  tags$hr(),
                                                  uiOutput("box_1_1")
                                                  
                                              ),
                                              
                                              fluidRow(
                                                  wellPanel(
                                                      h4("Dotplot for multiple gene DE visualization"),
                                                      
                                                      selectInput(inputId = "select_markers_dotplot", label = strong("Choose genes:"),choices = all_genes_ra, multiple = T, selected = fav_genes),
                                                      withSpinner(plotOutput("marker_dotplot")),
                                                      tags$hr(),
                                                      uiOutput("box_1_2")
                                                      
                                                      # )
                                                  )
                                              )
                                          )
                                          
                                      )
                                      
                             ),
                             tabPanel("Cell population view",
                                      wellPanel(
                                          
                                          fluidRow(
                                              column(
                                                  4,
                                                  uiOutput("cluster_ids"),
                                                  ##Comparing conditions new addition
                                                  selectInput(inputId = "ra_conds", label = strong("Choose conditions to compare:"),choices = conds)
                                              )
                                          ),
                                          fluidRow(
                                              
                                              column(6,
                                                     wellPanel(
                                                         h4("DE scatterplot"),
                                                         fluidRow(
                                                             withSpinner(plotOutput("cell_type_plot", click = clickOpts(id ="plot_click"))), dataTableOutput("click_info"),
                                                             
                                                             
                                                         ),
                                                         tags$hr(),
                                                         uiOutput("box_1_3a")
                                                         
                                                     )
                                              ),
                                              
                                              column(
                                                  6,
                                                  wellPanel(
                                                      h4("DE table"),
                                                      
                                                      # sliderInput(inputId = "top_genes", label = strong("Number of top DE genes:"), value = 100, min = 1, max = dim(cluster)[1], step = 1),
                                                      uiOutput("topdegenes"),
                                                      withSpinner(DTOutput("top_de_genes")),
                                                      tags$hr(),
                                                      uiOutput("box_1_3b")
                                                      
                                                      # )
                                                  )
                                              )
                                          )
                                          
                                      )
                                      
                             )
                         )
                )
                
)


##server function to compute the outputs
server = function(input, output) {
    
    
    ######################
    ### Selection of seurat object and plotting UMAP plots under subtitle 2.1
    ######################
    ##Selecting the Seurat object
    umap_clusters = reactive({
        alveri[[1]]
        # tcells_combined_umap_list_res[[1]]
    })
    
    ##Plotting UMAP plots for clustering
    output$all_groups = renderPlot({
        
        if(length(unique(umap_cluster_modified_ren_reo()[["FinalClusters"]][["FinalClusters"]])) == length(cluster.colours)){
            DimPlot(umap_cluster_modified_ren_reo(), reduction = "umap", split.by = "group", label.size = 6, order = T, cols = cluster.colours)
            
        } else {
            DimPlot(umap_cluster_modified_ren_reo(), reduction = "umap", split.by = "group", label.size = 6, order = T)
        }
        
    })
    
    
    ######################
    ### END Selection of seurat object and plotting UMAP plots under subtitle 2.1
    ######################
    
    
    ######################
    ### Cluster markers plots and tables under subtitle 2.2
    ######################
    ##Dynamic input field for selecting cluster to plot table of markers
    output$dyn_clusters <- renderUI({
        selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = unique(umap_clusters()[["FinalClusters"]][["FinalClusters"]]), multiple = F)
    })
    
    ##Displaying table of cluster markers for annotating cell types
    cluster_markers = reactive({
        # ra_macrophage_combined_clusters_tables_res[[1]][[(as.numeric(input$marker_genes_cluster) + 1)]] %>% rownames_to_column(var = 'genes') %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2))
        ra_macrophage_combined_clusters_tables_res[[1]][[(as.numeric(input$marker_genes_cluster) + 1)]] %>% rownames_to_column(var = 'genes') %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)) %>% dplyr::select(genes,contains(c("avg", "adj"))) %>% mutate_if(is.numeric, ~sprintf("%.3f", .))
        
        
    })
    
    output$topclgenes <- renderUI({
        sliderInput(inputId = "topclgenes_i", label = strong("Number of top cluster markers:"), value = 100, min = 1, max = dim(cluster_markers())[1], step = 1)
    })
    
    output$top_conserved_genes = DT::renderDataTable({
        #numeric_cols =  colnames(data.frame(cluster_markers()))[which_numeric_cols(data.frame(cluster_markers()))]
        
        #   # Javascript-enabled table.
        #   datatable(
        DT::datatable(
            head(cluster_markers(), n = input$topclgenes_i),
            selection = "single",
            rownames = FALSE,
            filter = list(position = "top", plain = TRUE),
            style = "default",
            extensions = c("Buttons","Scroller"),
            options = list(
                deferRender = TRUE,
                scrollY = 350,
                #scroller = TRUE,
                #lengthMenu = FALSE,
                autoWidth = FALSE,
                dom = "Blfrtip",
                buttons = 
                    list(list(
                        extend = "collection",
                        buttons = c("csv", "pdf"),
                        text = "Download"
                    )),  # end of buttons customization
                
                # customize the length menu
                lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), # end of lengthMenu customization
                pageLength = 10
            ), fillContainer = TRUE
        ) #%>%
        #DT::formatSignif(columns = numeric_cols, digits = 3)
    }, server = TRUE)
    
    
    ##Preparing and plotting UMAP cluster markers for annotating cell types
    umap_cluster_modified_rna = reactive({
        umap_cluster_modified_ul = umap_clusters()
        DefaultAssay(umap_cluster_modified_ul) = "RNA"
        umap_cluster_modified_ul
    })
    
    
    output$top_markers_umap <- renderUI({
        selectInput(inputId = "select_markers_umap", label = strong("Select marker to visualize in clusters:"), choices = cluster_markers()[,1], multiple = T, selected = head(cluster_markers()[,1], n=4))
        
    })
    output$conserved_markers_umap = renderPlot({
        FeaturePlot(umap_cluster_modified_rna(), features = input$select_markers_umap, min.cutoff = "q9")
    })
    
    
    
    ##information box
    output$box_2_2 =renderUI({
        wellPanel(style = "background:#385A4F",
                  tags$hr(),
                  tags$p(style = "font-family:Arial;color:white",
                         paste("Listing top cluster marker genes irrespective of cell state which can subsequently be used in labelling the cluster. Here, markers are genes highly expressed in a cluster as compared to all other clusters in both", conditions[1], conditions[2], conditions[3], conditions[4], "and", conditions[5],".")
                  ),
                  tags$hr()
        )
    })
    
    
    ######################
    ### END Cluster markers plots and tables under subtitle 2.2
    ######################
    
    
    ######################
    ### Labelling clusters under subtitle 2.3
    ######################
    ##Generating dynamic fields for labelling UMAP clusters and initializing the fields with placeholders
    output$cluster_annot <- renderUI({
        
        if(length(unique(umap_cluster_modified_rna()[["FinalClusters"]][["FinalClusters"]])) == length(cluster_names)){
            do.call(flowLayout, 
                    lapply(0:(length(unique(umap_cluster_modified_rna()[["FinalClusters"]][["FinalClusters"]]))-1), function(x) {
                        textInput(inputId = paste("labeller", x), label = strong(paste("Input cluster", x, "label")), value = cluster_names[x+1])
                    })
            )
            
        } else {
            do.call(flowLayout,
                    lapply(0:(length(unique(umap_cluster_modified_rna()[["FinalClusters"]][["FinalClusters"]]))-1), function(x) {
                        textInput(inputId = paste("labeller", x), label = strong(paste("Input cluster", x, "label")), value = x)
                    })
            )
        }
    })
    
    ##Storing names decided on by the researchers based on optimal clustering to start off the differential expression visualization
    annotation = reactiveValues(annot = cluster_names)
    
    ##Observer to allow updating of cluster names dynamically as typed in
    observe({
        
        req(unlist(lapply(0:(length(unique(umap_cluster_modified_rna()[["FinalClusters"]][["FinalClusters"]]))-1), function(x) {
            new_cluster_name = input[[paste("labeller",x)]]
        })))
        annotation$annot = unlist(lapply(0:(length(unique(umap_cluster_modified_rna()[["FinalClusters"]][["FinalClusters"]]))-1), function(x) {
            new_cluster_name = input[[paste("labeller",x)]]
        }))
    })
    
    ##Dynamic input for selecting celltypes (clusters) for diffential expression visualization
    output$cluster_ids <- renderUI({
        umap_names = annotation$annot
        #FindMarkers(stim_markers(), ident.1 = paste(input$select_cell_type, "KO", sep = "_"), ident.2 = paste(input$select_cell_type, "WT", sep = "_"), verbose = FALSE)
        
        if(length(unique(umap_cluster_modified_rna()[["FinalClusters"]][["FinalClusters"]])) == length(umap_names)){
            umap_names = annotation$annot
            selectInput(inputId = "select_cell_type", label = strong(paste("Select cell population to compare gene expression",":")),choices = umap_names, multiple = F)
            
        } else {
            selectInput(inputId = "select_cell_type", label = strong("Select cell population to compare gene expression across conditions:"),choices = unique(umap_cluster_modified_rna()[["FinalClusters"]][["FinalClusters"]]), multiple = F)
        }
        
    })
    
    ##Renaming clusters
    umap_cluster_modified_ren_reo = reactive({
        umap_names = annotation$annot
        umap_cluster_modified1 = umap_cluster_modified_rna()
        if(length(unique(umap_cluster_modified1$FinalClusters)) == length(umap_names)){
            names(umap_names) <- levels(umap_cluster_modified1)
            umap_cluster_modified1 <- RenameIdents(umap_cluster_modified1, umap_names)
            
        } else {
            umap_cluster_modified1
        }
        
    })
    
    # umap_cluster_modified_umap = reactive({
    #   DimPlot(umap_cluster_modified_ren_reo(), label = TRUE, label.size = 6)
    #   
    # })
    output$labelled_umap = renderPlot({
        
        if(length(unique(umap_cluster_modified_ren_reo()[["FinalClusters"]][["FinalClusters"]])) == length(cluster.colours)){
            DimPlot(umap_cluster_modified_ren_reo(), label = TRUE, label.size = 6, order = T, cols = cluster.colours)
            
        } else {
            DimPlot(umap_cluster_modified_ren_reo(), order = T)
        }
        
    })
    
    ######################
    ###END Labelling clusters under subtitle 2.3
    ######################
    
    
    ######################
    ### Differential expression dynamic using UMAP and Violin plots under subtitle 1.1
    ######################
    stim_markers = reactive({
        
        umap_cluster_modified = umap_cluster_modified_ren_reo()
        umap_cluster_modified$celltype.group <- interaction(Idents(umap_cluster_modified), umap_cluster_modified$group, sep = "_")
        umap_cluster_modified$celltype <- Idents(umap_cluster_modified)
        Idents(umap_cluster_modified) <- "celltype.group"
        umap_cluster_modified
    })
    
    
    #Functions to update differentially expressed genes
    output$de_stim_vs_ctrl_um = renderPlot({
        
        FeaturePlot(stim_markers(), features = input$de_genes, split.by = "group", max.cutoff = 3,cols = c("lightgrey", "blue"))
    })
    
    output$de_stim_vs_ctrl_vp = renderPlot({
        
        plots <- VlnPlot(stim_markers(), features = input$de_genes, split.by = "group", group.by = "celltype", pt.size = 0, combine = FALSE, multi.group = T, cols = group.cols, assay = "RNA") 
        for(i in 1:length(plots)) {
            plots[[i]] <- plots[[i]] + stat_summary(fun.y= median, geom='point', size = 2, colour = "black", position = position_dodge(0.9))
        }
        CombinePlots(plots)
        
    })
    
    ## Define two panels for UMAP and violin plots
    
    
    ##Information box
    output$box_1_1 <- renderUI({
        wellPanel(style = "background:#385A4F",
                  tags$hr(),
                  tags$p(style = "font-family:Arial;color:white",
                         #titlePanel("includeText"),
                         #fluidRow(
                         #column(12,
                         paste("Comparison of", input$de_genes, "expression between using violin plots and umap.")
                         #))
                  ),
                  tags$hr()
        )  })
    
    ######################
    ### END Differential expression dynamic using UMAP and Violin plots under subtitle 1.1
    ######################
    
    
    ######################
    ### Differential expression using dotplot under subtitle 1.2
    ######################
    ##Dotplot for DE comparison between KO and WT across cell types
    output$marker_dotplot = renderPlot({
        umap_cluster_modified_ren_reo = umap_cluster_modified_ren_reo()
        umap_cluster_modified_ren_reo@meta.data$grp_od <- umap_cluster_modified_ren_reo@meta.data$group
        umap_cluster_modified_ren_reo@meta.data <- umap_cluster_modified_ren_reo@meta.data %>% mutate(grp_od = case_when(grp_od == "Healthy" ~ 4,grp_od == "UPA" ~ 3,grp_od == "Naive RA" ~ 2,grp_od == "Resistant RA" ~ 1,grp_od == "Remission RA" ~ 0))
        
        ##Now you can arrange the data in increasing or decreasing order of genotype using the GEN column to get the dotplot the way you want.
        umap_cluster_modified_ren_reo@meta.data <- dplyr::arrange(umap_cluster_modified_ren_reo@meta.data, umap_cluster_modified_ren_reo@meta.data$grp_od)
        
        ##or **(not the '-' sign)**
        umap_cluster_modified_ren_reo@meta.data <- dplyr::arrange(umap_cluster_modified_ren_reo@meta.data, -umap_cluster_modified_ren_reo@meta.data$grp_od)
        DotPlot(umap_cluster_modified_ren_reo, features = input$select_markers_dotplot, cols = group.cols, dot.scale = 6, split.by = "group") + RotatedAxis()  + coord_flip()
    })
    
    ##Information box
    output$box_1_2 <- renderUI({
        wellPanel(style = "background:#385A4F",
                  tags$hr(),
                  tags$p(style = "font-family:Arial;color:white",
                         
                         paste("Comparison of gene expression between", conditions[1], conditions[2], conditions[3], conditions[4], "and", conditions[5], "cells across clusters using a dotplot. The genes are on the y-axis and the clusters on the x-axis. Strong-green, Orange, Dark-orange1, Red and Blue representing healthy, UPA, naive RA, resistant RA and remission RA respectively with the increase in intensity of the respective colour correlating with the increase in the average level of gene expression across all cells in the cluster. The size of the dot corresponds to the percentage of cells in the cluster expressing the gene.")
                         
                  ),
                  tags$hr()
        )})
    ######################
    ###END Differential expression using dotplot under subtitle 1.2
    ######################
    
    
    ######################
    ### Differential expression using ggplot and tables under subtitle 1.3
    ######################
    
    ##Retrieving table for DE expression from precomputed list
    genes_in_de_order = reactive({
        umap_names = annotation$annot
        
        if(length(unique(umap_cluster_modified_rna()[["FinalClusters"]][["FinalClusters"]])) == length(umap_names)){
            ra_macrophage_combined_de_tables[[1]][[as.numeric(match(input$select_cell_type,umap_names))]][[as.numeric(match(input$ra_conds,conds))]]
            
        } else {
            # ra_macrophage_combined_de_tables[[1]][[(as.numeric(input$select_cell_type) + 1)]]
            ra_macrophage_combined_de_tables[[1]][[(as.numeric(input$select_cell_type) + 1)]][[as.numeric(match(input$ra_conds,conds))]]
            
        }
        
        
        
    })
    
    ##Retrieving table for DE expression from precomputed list
    top_de_g = reactive({
        genes_in_de_order()  %>% rownames_to_column(var = 'genes') %>% filter(p_val_adj <= 0.05) %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)) %>% mutate_if(is.numeric, ~sprintf("%.3f", .))
    })
    
    
    output$topdegenes <- renderUI({
        sliderInput(inputId = "topdegenes_i", label = strong("Number of top DE genes:"), value = 100, min = 1, max = dim(top_de_g())[1], step = 1)
    })
    
    ##for table
    # top_de_g_tbl = reactive({
    #   head(genes_in_de_order(), n = input$top_de_genes)
    #   })
    output$top_de_genes = DT::renderDataTable({
        # numeric_cols =
        #   colnames(data.frame(top_de_g()))[which_numeric_cols(data.frame(top_de_g()))]
        
        #   # Javascript-enabled table.
        #   datatable(
        DT::datatable(
            head(top_de_g(), n = input$topdegenes_i),
            # colnames = c(
            #   "Gene",
            #   "scRNA-seq cluster",
            #   "Avg_logFC",
            #   "% in interest cluster",
            #   "% in other clusters",
            #   "Adjusted P"
            # ),
            selection = "single",
            rownames = FALSE,
            filter = list(position = "top", plain = TRUE),
            style = "default",
            extensions = c("Buttons","Scroller"),
            options = list(
                deferRender = TRUE,
                scrollY = 350,
                #scroller = TRUE,
                #lengthMenu = FALSE,
                autoWidth = FALSE,
                dom = "Blfrtip",
                buttons = 
                    list(list(
                        extend = "collection",
                        buttons = c("csv", "pdf"),
                        text = "Download"
                    )),  # end of buttons customization
                
                # customize the length menu
                lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), # end of lengthMenu customization
                pageLength = 10
            )
        ) #%>%
        # DT::formatSignif(columns = numeric_cols, digits = 3)
    }, server = TRUE)
    
    ##Allowing for download of DE table
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("differentially_expressed_genes_in",input$select_cell_type,input$ra_conds, ".csv", sep = "_")
            
        },
        content = function(file) {
            write.csv(genes_in_de_order() %>% rownames_to_column(var = 'genes') %>% filter(p_val_adj <= 0.05)  %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)), file)
        }
    )
    
    ##Retrieving table for DE scatterplotfrom precomputed list
    cell_type_de = reactive({
        umap_names = annotation$annot
        
        if(length(unique(umap_cluster_modified_rna()[["FinalClusters"]][["FinalClusters"]])) == length(umap_names)){
            ra_macrophage_combined_de_ggplots_table[[1]][[as.numeric(match(input$select_cell_type,umap_names))]]
            
        } else {
            ra_macrophage_combined_de_ggplots_table[[1]][[(as.numeric(input$select_cell_type) + 1)]]
        }
        
    })
    
    ##preparing ggplot for average DE expression for genes above
    cell_type_de_plot = reactive({
        #theme_set(theme_cowplot())
        grob <- grobTree(textGrob("Click on points to diplay more information about the gene", x=0.1,  y=0.95, hjust=0,
                                  gp=gpar(col="red", fontsize=13, fontface="italic")))
        ggplot(data=cell_type_de(), aes_string(paste("`",unlist(str_split(input$ra_conds[1], " VS "))[1],"`", sep=""), paste("`",unlist(str_split(input$ra_conds[1], " VS "))[2],"`", sep=""))) + geom_point() + ggtitle(input$select_cell_type) + annotation_custom(grob) + theme_bw()
        
        
    })
    
    ##plotting ggplot for average DE expression for genes above
    output$cell_type_plot = renderPlot({
        #theme_set(theme_cowplot())
        cell_type_de_plot()
    })
    
    
    ##Displaying further details upon clicking points
    displayed_text <- reactive({
        req(input$plot_click)
        nearPoints(cell_type_de(), input$plot_click)
        
    })
    
    ##Displaying table with gene details upon click of point in DE scatterplot
    output$click_info <- renderDataTable({
        req(displayed_text())
        
        DT::datatable(displayed_text(),
                      extensions=c('Scroller'),
                      options = list(dom = 'Bfrtip',
                                     scroller = TRUE,
                                     scrollX=TRUE)) 
    }, escape = F)
    
    ##Information box
    output$box_1_3a <- renderUI({
        wellPanel(style = "background:#385A4F",
                  tags$hr(),
                  tags$p(style = "font-family:Arial;color:white",
                         paste("Comparison of average gene expression between", input$ra_conds, input$select_cell_type, "using a scatter plot.")
                  ),
                  tags$hr()
        )
    })
    
    output$box_1_3b <- renderUI({
        
        wellPanel(style = "background:#385A4F",
                  tags$hr(),
                  tags$p(style = "font-family:Arial;color:white",
                         paste("Listing differentially expressed genes (adjusted P value <0.05) in" , input$ra_conds, "among", input$select_cell_type,".")
                  ),
                  tags$hr()
        )
    })
    ######################
    ### END Differential expression using ggplot and tables under subtitle 1.3
    ######################
    
    
}


shinyApp(ui = ui, server = server)


