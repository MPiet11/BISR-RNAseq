VolcanoplotFunc <- function(toptable, NominalCutoff, AdjustedCutoff, LabellingCutoff, FCCutoff, main, type)
{
  toptable$Significance <- "NS"
  toptable$Significance[(abs(toptable$logFC) > FCCutoff)] <- "FC"
  toptable$Significance[(toptable$adj.P.Val<AdjustedCutoff)] <- "FDR"
  toptable$Significance[(toptable$adj.P.Val<AdjustedCutoff) & (abs(toptable$logFC)>FCCutoff)] <- "FC_FDR"
  table(toptable$Significance)
  toptable$Significance <- factor(toptable$Significance, levels=c("NS", "FC", "FDR", "FC_FDR"))


  plot <- ggplot(toptable, aes(x=logFC, y=-log10(adj.P.Val))) +

    #Add points:
    #   Colour based on factors set a few lines up
    #   'alpha' provides gradual shading of colour
    #   Set size of points
    geom_point(aes(color=factor(Significance)), alpha=1/2, size=0.8) +

    #Choose which colours to use; otherwise, ggplot2 choose automatically (order depends on how factors are ordered in toptable$Significance)
    scale_color_manual(
      values = c(
        NS = "grey30",
        FC = "forestgreen",
        FDR = "royalblue",
        FC_FDR = "red2"
      ),
      labels = c(
        NS = "Not Significant",
        FC = paste("LogFC>|", FCCutoff, "|", sep = ""),
        FDR = paste("Adj.P.val <", AdjustedCutoff, sep = ""),
        FC_FDR = paste("Adj.P.val <", AdjustedCutoff, " & LogFC > ", FCCutoff, " ", sep =
                         "")
      )
    ) +

    #Set the size of the plotting window
    theme_bw(base_size=24)  +

    #Modify various aspects of the plot text and legend
    theme(legend.background=element_rect(),
          plot.title=element_text(angle=0, size=12, face="bold", vjust=1),

          panel.grid.major=element_blank(),   #Remove gridlines
          panel.grid.minor=element_blank(),   #Remove gridlines

          axis.text.x=element_text(angle=0, size=12, vjust=1),
          axis.text.y=element_text(angle=0, size=12, vjust=1),
          axis.title=element_text(size=12),

          #Legend
          legend.position="top",          #Moves the legend to the top of the plot
          legend.key=element_blank(),     #removes the border
          legend.key.size=unit(0.5, "cm"),    #Sets overall area/size of the legend
          legend.text=element_text(size=8),   #Text size
          title=element_text(size=8),     #Title text size
          legend.title=element_blank()) +     #Remove the title

    #Change the size of the icons/symbols in the legend
   guides(colour = guide_legend(override.aes=list(size=2.5))) +

   #  #Set x- and y-axes labels

   #xlab(label_bquote(~-log[10]~adjusted~italic(P)))
   #ylab(label_bquote( ~ log[2] ~ "fold change"))
   xlab("log2 fold change") +
   ylab("-log10 adjusted pvalue") +

  # if (tmp.data[1]!="no value") {
  #   p <- p + geom_point() + geom_line()
  # } else {
  #   p <- p + geom_line()
  # }



    #Set the axis limits
    #xlim(-6.5, 6.5) +
    #ylim(0, 100) +

    #Set title
   ggtitle(main) +

    #Tidy the text labels for a subset of genes
    # geom_text(data=subset(toptable, FDR<LabellingCutoff & abs(logFC)>FCCutoff),
    #           aes(label=rownames(subset(toptable, FDR<LabellingCutoff & abs(logFC)>FCCutoff))),
    #           size=2.25,
    #           #segment.color="black", #This and the next parameter spread out the labels and join them to their points by a line
    #           #segment.size=0.01,
    #           check_overlap=TRUE,
    #           vjust=1.0) +

    #Add a vertical line for fold change cut-offs
    geom_vline(xintercept=c(-FCCutoff, FCCutoff), linetype="longdash", colour="black", size=0.4) +

    #Add a horizontal line for P-value cut-off
    geom_hline(yintercept=-log10(AdjustedCutoff), linetype="longdash", colour="black", size=0.4)

  #return(plot)

  # return(ggplotly(plot))

  if (type == "static") {
    return(plot)
  }
  if(type == "interactive") {
    return(ggplotly(plot))
  }


}
results <- topTable(fit2, coef = 1, n = Inf)
CZ10FF_vs_1E_toptable <- results %>% rownames_to_column("gene_name")%>%
  select(gene_name, logFC, adj.P.Val)
VolcanoplotFunc(
  toptable = CZ10FF_vs_1E_toptable,
  AdjustedCutoff = 0.05,
  LabellingCutoff = 0.05,
  FCCutoff = 2.0,
  main = "EdgeR results",
  type = "static"
)

#
#
# dat <- data.frame(cond = rep(c("A", "B"), each=10),
#                   xvar = 1:20 + rnorm(20,sd=3),
#                   yvar = 1:20 + rnorm(20,sd=3))
# is(dat)
# p <- ggplot(dat, aes(x=xvar, y=yvar)) +
#   geom_point(shape=1)      # Use hollow circles
#
# ggplotly(p)
#
# #fit2
# results <- topTable(fit2, coef = 1, n = Inf)
# CZ10FF_vs_1E_toptable <- results %>% rownames_to_column("gene_name")%>%
#   select(gene_name, logFC, adj.P.Val)
# # This works
#
#
# # Rounding values messing up Volcano plot
# CZ10FF_vs_1E_rounded <- as.tibble(lapply(CZ10FF_vs_1E_toptable, function(y)
#   ifelse(
#     is.na(as.numeric(y)), y, round(as.numeric(y), 2)
#   )))
# # rounded values mess the volcano plot
# EnhancedVolcanoEdgeR(
#   toptable = CZ10FF_vs_1E_rounded,
#   AdjustedCutoff = 0.05,
#   LabellingCutoff = 0.05,
#   FCCutoff = 2.0,
#   main = "EdgeR results"
# )
#
#
# head(CZ10FF_vs_1E_toptable)
# head(datalist$CZ10FF_vs_CZ1E)
# #colnames(results) <- c("logFC","AveExpr", "t", "P.Value", "FDR","B")
# #colnames(results)
# results <- rshiny_input$volcanto_plot %>%
#   select(gene_name, CZ10FF_vs_CZ1E_logFC, CZ10FF_vs_CZ1E_adj.P.Val)
# colnames(results) <- c("gene_name","logFC", "adj.P.Val")
#
# EnhancedVolcanoEdgeR(
#   toptable = results,
#   AdjustedCutoff = 0.05,
#   LabellingCutoff = 0.05,
#   FCCutoff = 2.0,
#   main = "EdgeR results"
# )
#
# ggplotly(
# ggplotly(p = ggplot2::last_plot(), width = NULL, height = NULL,
#          tooltip = "all", dynamicTicks = FALSE, layerData = 1,
#          originalData = TRUE, source = "A" )
#
#
# ggiris <- qplot(Petal.Width, Sepal.Length, data = iris, color = Species)
# ggplotly(ggiris)
#
#
# p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() +
#   facet_grid(vs ~ ., labeller = label_bquote(alpha ^ .(vs))) +
#   facet_grid(. ~ vs, labeller = label_bquote(cols = .(vs) ^ .(vs))) +
#   facet_grid(. ~ vs + am, labeller = label_bquote(cols = .(am) ^ .(vs)))
#
# ggplotly(p)
#
# # rshiny_input$volcanto_plot$CZ10FF_vs_CZ1E_logFC
# # rshiny_input$volcanto_plot$CZ10FF_vs_CZ1E_adj.P.Val
# # rshiny_input$volcanto_plot$gene_name
#
#
#
