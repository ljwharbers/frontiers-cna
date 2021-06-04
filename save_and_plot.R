require("ggplot2")
require("gridExtra")
require("cowplot")
require("ggpubr")
require("ggthemes")
require("ggsci")
require("viridis")

theme_set(theme_cowplot())

# Function to save a plot in different file formats and show it within the notebook
save_and_plot <- function(x, bname, width, height, dpi=300, 
                          use.cowplot=FALSE, ncol=1, nrow=1, base_height=4, base_width=NULL, 
                          base_aspect_ratio=1, plot=FALSE, onefile = TRUE){
  if( !use.cowplot ){
    png(filename=file.path(paste0(bname, ".png")), units="in", 
        height=height, width=width, res=dpi)
    print(x)
    while(length(dev.list())>0) invisible(dev.off())
    cairo_ps(filename=file.path(paste0(bname, "_cairo_ps.eps")), 
             onefile = TRUE, height=height, width=width, family="Helvetica", 
             pointsize=8, antialias="none")
    print(x)
    while(length(dev.list())>0) invisible(dev.off())
    cairo_pdf(filename=file.path(paste0(bname, "_cairo_pdf.pdf")), 
              onefile = TRUE, height=height, width=width, family="Helvetica", 
              pointsize=8, antialias="none")
    print(x)
    while(length(dev.list())>0) invisible(dev.off())
    postscript(file=file.path(paste0(bname, "_postscript.eps")), 
               onefile = TRUE, paper="special", height=height, width=width, 
               family="Helvetica", pointsize=8, horizontal=FALSE)
    print(x)
    # while(length(dev.list())>0) invisible(dev.off())
    # svg(file=file.path(paste0(bname, ".svg")), 
    #            onefile = TRUE, height=height, width=width, 
    #            family="Helvetica", pointsize=8)
    # print(x)
    # while(length(dev.list())>0) invisible(dev.off())
  }else{
    save_plot(x, filename=file.path(paste0(bname, ".png")), ncol = 
                ncol, nrow = nrow, base_height = base_height, base_width = base_width, 
              base_aspect_ratio = base_aspect_ratio, dpi=dpi)
    while(length(dev.list())>0) invisible(dev.off())
    save_plot(x, filename=file.path(paste0(bname, "_cairo_ps.eps")), 
              ncol = ncol, nrow = nrow, base_height = base_height, base_width = 
                base_width, base_aspect_ratio = base_aspect_ratio, device=cairo_ps)
    save_plot(x, filename=file.path(paste0(bname, "_cairo_pdf.pdf")), 
              ncol = ncol, nrow = nrow, base_height = base_height, base_width = 
                base_width, base_aspect_ratio = base_aspect_ratio, device=cairo_pdf)
    while(length(dev.list())>0) invisible(dev.off())
    save_plot(x, filename=file.path(paste0(bname, "_postscript.eps")), 
              ncol = ncol, nrow = nrow, base_height = base_height, base_width = 
                base_width, base_aspect_ratio = base_aspect_ratio, device="ps")
    while(length(dev.list())>0) invisible(dev.off())
    # save_plot(x, filename=file.path(paste0(bname, ".svg")), 
    #           ncol = ncol, nrow = nrow, base_height = base_height, base_width = 
    #             base_width, base_aspect_ratio = base_aspect_ratio, device="svg")
    # while(length(dev.list())>0) invisible(dev.off())
  }
  if( plot )
    print(x)
  dev.off()
}