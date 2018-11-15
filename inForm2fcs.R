#
# inForm2fcs.R
#
# Utilities for working with data tables produced by inForm software (PerkinElmer)
# for cell segmentation, phenotyping, quantification, and localization of Vectra scans.
# *_cell_seg_data.txt
# *_cell_seg_data_summary.txt
#
# 
# make_fcs()
# Enable flow-style analysis of Vectra tissue scans for inForm segmentation results.
# Read data produced by inForm software,
# write out fcs files, grouped by "Tissue Category" and "Phenotype" per each "Slide ID".
# Also want an "all" file per image. Including all phenotypes and tissue categories.
#
# "Sample Name" = "Slide ID" + x,y coords of the box drawn on the image.
#
# addZeros() - If there are 0 cells for a category for all regions in a region, the category is omitted. Add it with a count of 0.
# addZeros_inform_merged(fname)  make sure fname is not overridden
#
#
#  devtools::install_github("PerkinElmer/phenoptr", build_vignettes=TRUE); # not necessary
#

# start shiny app:  runApp("inFormA.R")

#library("phenoptr"); # for read_cell_seg_data(fname) # not necessary # just as good: csd <- read.table(fname, header=TRUE, sep="\t", as.is=TRUE, comment.char="", check.names=FALSE);
library("flowCore"); # for write.FCS()
library(shiny); # for UI

test_function <- function(arg1){
   print("test_function: %s", arg1);
}

#fname <- "DATA/120517 Panel 1 Batch Merge_cell_seg_data.txt"; outdir <- ".";
#make_fcs <- function(fname = "DATA/MTC Panel 1 021618 Batch Merge_cell_seg_data.txt", outdir=".", dropNormal=FALSE){
#make_fcs <- function(fname = "DATA/MTC Panel 4 (15-2) Thyroid Batch merge_cell_seg_data.txt", outdir=".", dropNormal=FALSE){
#make_fcs <- function(fname = "DATA/MTC Panel 4 (15-2) TILN Batch Merge_cell_seg_data.txt", outdir=".", dropNormal=FALSE){
make_fcs <- function(fname = "concat5_cell_seg_data.txt", outdir=".", dropNormal=FALSE, origfname=c()){
   t0 <- Sys.time();
   if(is.null(origfname)){
      origfname <- fname;
   }
   
   ####  ignore outdir  ###
   # Outfile names will start with outdir, 
   # then basename of infile,
   # convert spaces to underscores,
   # drop cell_seg_data.txt
   #outbase <- sprintf("%s/%s", outdir, 
   #   sub(pattern="_cell_seg_data.txt$", replacement="_", x=
   #      gsub(pattern=" ", replacement="_", x=basename(fname), fixed=TRUE)));
   outbase <- sub(pattern="_cell_seg_data.txt$", replacement="_", x=
         gsub(pattern=" ", replacement="_", x=basename(origfname), fixed=TRUE));

   print(sprintf("fname:%s", fname), q=F);
   #print(sprintf("outbase:%s", outbase), q=F);

   #csd <- read_cell_seg_data(fname); #  class(csd):   "tbl_df"     "tbl"        "data.frame"
   # just as good 
   csd <- read.table(fname, header=TRUE, sep="\t", as.is=TRUE, comment.char="", check.names=FALSE);
   colnames(csd) <- gsub(pattern=" (Normalized Counts, Total Weighting)", replacement="", x=colnames(csd), fixed=TRUE);


   if(dropNormal){
      # Drop rows if the "Path" column (2) contains the string "normal" or "Normal"
      #grepNormal_i <- grep("normal", csd$Path, ignore.case=TRUE, value=FALSE, fixed=TRUE);
      grepNormal_i <- grep("normal", csd$Path, ignore.case=TRUE, value=FALSE, fixed=FALSE);
      if(length(grepNormal_i) > 0){
         print(sprintf("Dropping %i normals from %s", length(grepNormal_i), basename(origfname)), q=F);
         print(unique(csd$Path[grepNormal_i]), q=F);
         csd <- csd[-grepNormal_i,];
      }
   }

   if(! ("Tissue Category" %in% colnames(csd)) ){
      csd <- cbind(csd, "Tissue Category"="single_tissue_category");
   }
   # How many cells of each phenotype are in each tissue category?
   print(table(csd$`Tissue Category`, csd$Phenotype), q=F);

   csd_mat <- as.matrix(csd);
   #  [1] "tag"      NOTE: this column is only sometimes included; can't rely on column numbers
   #  [2] "Path"                                        
   #  [3] "Sample Name"     :  for each box; n=294, ex: "111317 S09-8015 D Panel 12-2_[56051,22921].im3"; name inculdes image_filename and x,y coords
   #  [4] "Tissue Category" :  "Non-tumor" "Tumor", possibly "other"?                            
   #  [5] "Phenotype"       :  "Other" "CK+"   "CD4"   "Treg"  "CD8"                             
   #  [6] "Cell ID"                   
   
   #  unique combinations of "Sample Name" and "Tissue Category": 586
   #  unique combinations of "Sample Name", "Tissue Category", "Phenotype": 2105

   ####   subsample for testing   ####
   #samplesize <- 10000;
   #samplesize <- 100000;
   #samplesize <- 1000;
   #ris <- sample.int(n=nrow(csd), size=samplesize);
   #csd_mat <- csd_mat[ris,];

   # Drop columns we don't need.
   allMeanCols <- grep("Mean", colnames(csd_mat), value=TRUE, fixed=TRUE);
   colsNeeded <- c(allMeanCols, "Slide ID", "Phenotype", "Tissue Category"); 
   wcols <- which(colnames(csd_mat) %in% colsNeeded);
   csd_mat <- csd_mat[,wcols];
   # csd_mat is character, has columns: allMeanCols, "Slide ID", "Phenotype", "Tissue Category"

   # Which are mean columns? These will be converted from char to numeric and written to fcs out.
   wmc <- which(colnames(csd_mat) %in% allMeanCols);

   #gotOne <- FALSE; # testing...

   # for each image (Slide ID) 
   for(slide_img_name in unique(csd_mat[,"Slide ID"])){

      #if(gotOne){
      #   break; # testing...
      #}

      slide_out_base <- sprintf("%s%s", outbase, gsub(pattern=" ", replacement="_", x=slide_img_name, fixed=TRUE));
      wslide <- which(csd_mat[,"Slide ID"] == slide_img_name);

      # Make the "all" file. All phenotypes and tissue types for this image.
      doublematrix <- matrix(as.double(csd_mat[wslide, wmc]), nrow=length(wslide), byrow=FALSE);
      colnames(doublematrix) <- colnames(csd_mat)[wmc];
      ff <- flowFrame(doublematrix);
      img_out_filename <- sprintf("%s_all.fcs", slide_out_base);
      print(sprintf("%s  %i cells", img_out_filename, length(wslide)), q=F);
      write.FCS(ff, filename=img_out_filename);

      for(tissueCat in unique(csd_mat[,"Tissue Category"])){
         wtc <- which(csd_mat[,"Tissue Category"] == tissueCat);

         # all phenotypes, by tissue category
         w_slide_t <- intersect(wslide, wtc); # which slide and tissue category
         if(length(w_slide_t) > 0){
            doublematrix <- matrix(as.double(csd_mat[w_slide_t, wmc]), nrow=length(w_slide_t), byrow=FALSE);
            colnames(doublematrix) <- colnames(csd_mat)[wmc];
            ff <- flowFrame(doublematrix);
            out_filename <- sprintf("%s_%s_all.fcs", slide_out_base, tissueCat);
            print(sprintf("%s  %i cells", out_filename, length(w_slide_t)), q=F);
            write.FCS(ff, filename=out_filename);
         } else{
            print(sprintf("skipping %s all %s  0 cells", slide_out_base, tissueCat), q=F);
         }

         for(pheno in unique(csd_mat[,"Phenotype"])){
            wph <- which(csd_mat[,"Phenotype"] == pheno);
         
            w_slide_ph_t <- intersect(wslide, intersect(wph, wtc));
            if(length(w_slide_ph_t) > 0){
               doublematrix <- matrix(as.double(csd_mat[w_slide_ph_t, wmc]), nrow=length(w_slide_ph_t), byrow=FALSE);
               colnames(doublematrix) <- colnames(csd_mat)[wmc];
               ff <- flowFrame(doublematrix);
               #out_filename <- sprintf("%s_%s_%s.fcs", slide_out_base, pheno, tissueCat);
               # Drop / from file name. write.FCS can't handle it.
               out_filename <- sprintf("%s_%s_%s.fcs", slide_out_base, gsub("/","-",pheno), tissueCat);
               print(sprintf("%s  %i cells", out_filename, length(w_slide_ph_t)), q=F);
               write.FCS(ff, filename=out_filename);
            } else{
               print(sprintf("skipping %s %s %s  0 cells", slide_out_base, pheno, tissueCat), q=F);
            }

         } # for each Tissue Category
      } # for each Phenotype 
   } # for each Slide ID (image)
   t1 <- Sys.time();
   print(t1-t0);
} # make_fcs()


# Manual quality control.
test_readback <- function(){
   fcsname <- "120517_Panel_1_Batch_Merge_111317_09-MS-9102-26_Panel_12-2_Other_Non-tumor.fcs"
   fcsname <- "120517_Panel_1_Batch_Merge_111317_JS13-5443_D1_Panel_12-2_CD4_Non-tumor.fcs"
   rb <- read.FCS(fcsname, transformation=c());


   fcsname <- "csd/5.9.18_Merged_cell_Segmentation_Data_1_CD8+ T cell_Lamina Propria.fcs";
   rd <- read.FCS(fcsname, transformation=c());
   fcsname <- "css/5.9.18_Merged_cell_Segmentation_Data_1_CD8+ T cell_Lamina Propria.fcs";
   rs <- read.FCS(fcsname, transformation=c());
}


testing <- function(){

   #hist(csd$`Entire Cell PD-1 (Opal 570) Mean`, breaks=400);
   #hist(csd$`Entire Cell PD-1 (Opal 570) Mean`[ris], breaks=400);
   #hist(log(csd$`Entire Cell PD-1 (Opal 570) Mean`), breaks=400);
   hist(log(csd$`Entire Cell PD-1 (Opal 570) Mean`[ris]), breaks=1000);

   hist(log(csd$`Entire Cell PD-1 (Opal 570) Mean`), breaks=10000);
   hist((csd$`Entire Cell PD-1 (Opal 570) Mean`), breaks=10000);
   hist((csd$`Entire Cell PD-1 (Opal 570) Mean`), breaks=10000, xlim=c(-.1,1));

   #csd <- csd %>% mutate(pdl1_plus=`Entire Cell PDL1 (Opal 520) Mean`>3)
   #csd <- csd %>% mutate(pdl1_plus=`Entire Cell PD-1 (Opal 570) Mean` > 1);
   #table(csd$pdl1_plus, csd$Phenotype)

   allMeanCols <- grep("Mean", colnames(y), value=TRUE, fixed=TRUE);
   nucMeanCols <- grep("Nucleus.*Mean", colnames(y), value=TRUE, fixed=FALSE);
   MembraneMeanCols <- grep("Membrane.*Mean", colnames(y), value=TRUE, fixed=FALSE);
   wholeCellMeanCols <- grep("Entire Cell.*Mean", colnames(y), value=TRUE, fixed=FALSE);

   # split based on phenotype or tissue type?
   # "Tisue Category", "Sample Name", "Phenotype"

   wmc <- which(colnames(y) %in% allMeanCols);
   doublematrix <- matrix(as.double(y[,wmc]), nrow=nrow(y), byrow=FALSE);
   #doublematrix <- matrix(as.double(y[1:5,wmc]), nrow=5, byrow=FALSE);
   colnames(doublematrix) <- colnames(y)[wmc];
   fmat <- flowFrame(as.numeric(y[, wmc]));

   summary_fname <- "DATA/120517 Panel 1 Batch Merge_cell_seg_data_summary.txt";

   dev.new();
   plot(as.numeric(csd$`Entire Cell PD-1 (Opal 570) Mean`[ris]), as.numeric(csd$`Membrane CD8 (Opal 650) Mean`[ris]), pch=".");
   plot(log(as.numeric(csd$`Entire Cell PD-1 (Opal 570) Mean`[ris])), log(as.numeric(csd$`Membrane CD8 (Opal 650) Mean`[ris])), pch=".");
   plot(log(as.numeric(csd$`Entire Cell PD-1 (Opal 570) Mean`[ris])), log(as.numeric(csd$`Membrane CD4 (Opal 620) Mean`[ris])), pch=".");
   pf <- levels(as.factor(csd$Phenotype[ris]));
   pfl <- levels(as.factor(csd$Phenotype[ris]));
   pfcols <- rep(1,length(csd$Phenotype[ris]));
   plot(log(as.numeric(csd$`Membrane CD8 (Opal 650) Mean`[ris])), log(as.numeric(csd$`Membrane CD4 (Opal 620) Mean`[ris])), pch=".", col=as.factor(csd$Phenotype[ris]));
   # add color based on phenotype

}


testing_addZeros <- function(fname = "DATA/121217 Panel 1 Repeat Merge_cell_seg_data_summary.txt", outdir="."){
   csd <- read_cell_seg_data(fname); #  class(csd):   "tbl_df"     "tbl"        "data.frame"
   csd_mat <- as.matrix(csd);
   usns <- unique(csd$'Sample Name');
   counts <- rep(0, length(usns));
   si <- 0;
   for(sn in unique(csd$'Sample Name')){
      if(length(wsn) != 18){
         print(sn);
      }
      si <- si + 1;
      wsn <- which(csd$'Sample Name' == sn);
      counts[si] <- length(wsn);
   }
   
   # copy the line for Tissue Category and Phenotype == All
   # fill numerics with 0
}


# Split output from addZeros() into seperate file per "Patient Tissue ID"
# called by addZeros()
cell_seg_summary_per_patient <- function(new_css){
   all_ids <- unique(new_css$'Sample Name');
   for(thisid in all_ids){
      wi <- which(new_css$'Sample Name' == thisid);
      subcss <- new_css[wi,];
      outfile <- sprintf("%s_cell_seg_data_summary.txt", thisid);
      write.table(subcss, file=outfile, quote=FALSE, sep="\t", row.names=FALSE); 
   }
}


# Inform software omits rows for Phenotype categories with 0 cells in all samples.
# Add entries for zeros so all samples have an entry for each Phenotype for each Tissue Category.
# fname may be from upload to tmp directory. If origfname is supplied, put output with that.
addZeros_inform_merged <- function(fname=c(), origfname=c()){
   # example
   #fname <- "/Volumes/HI3-Microscope/Data/Vectra3/Diamond Lab/Panel 41 HLADR Export/Panel 41 MERGE file_cell_seg_data_summary.txt";

   #print(fname, q=F);

   css <- read.table(fname, header=TRUE, sep="\t", as.is=TRUE, comment.char="", check.names=FALSE);
   added_TC <- FALSE;
   if(!('Tissue Category' %in% colnames(css))){
      css <-  cbind(css, 'Tissue Category'="All")
      added_TC <- TRUE;
   }
   # Check for and remove extra header lines left from sloppy manual mergers.
   wheaderlines <- which(css$'Sample Name' == "Sample Name");
   if(length(wheaderlines) > 0){
      css <- css[-wheaderlines,];
   }

   new_css <- css;
   all_cat <- unique(css$'Tissue Category');
   all_ph <- unique(css$Phenotype);
   all_samplenames <- unique(css$'Sample Name');
   print(sprintf("%i Tissue Categories:", length(all_cat)), q=F);
   print(paste(all_cat, collapse="|"), q=F);
   print(sprintf("%i Phenotypes:", length(all_ph)), q=F);
   print(paste(all_ph, collapse="|"), q=F);
   print(sprintf("%i Sample Names", length(all_samplenames)), q=F);
   print(sprintf("Expect %i events.", length(all_cat)*length(all_ph)*length(all_samplenames)), q=F);
   for(this_sample in all_samplenames){
      #wis <- which(this_sample == css[,3]);
      wis <- which(this_sample == css$'Sample Name');
      #print(length(wis));
      thiscss <- css[wis,];
      alli <- which(thiscss$Phenotype == "All" & thiscss$'Tissue Category' == "All");
      if(length(alli) != 1){
         #print(thiscss);
         print(this_sample);
         print(length(wis));
         print(alli);
         print(length(alli));
         stop("problem");
      }
      if(!(all(all_cat %in% thiscss$'Tissue Category'))){
         print(thiscss);
         stop(" didn't find all tissue categories");
      }
      # For each tissue category
      for(tc in all_cat){
         alli_tc <- which(thiscss$Phenotype == "All" & thiscss$'Tissue Category' == tc);
         if(length(alli_tc) != 1){
            print(tc);
            print(alli_tc);
            print(length(alli_tc));
            print(thiscss[,c(1,4:8)]);
            break;
            stop(" didn't find all for this category");
         }
         wtc <- which(thiscss$'Tissue Category' == tc);
         thiscss_tc <- thiscss[wtc,]; # subset of rows for this Phenotype and Tissue Category
         # For each phenotype
         for(ph in all_ph){ 
            if(!(ph %in% thiscss_tc$Phenotype)){
               # add it. Start with "All" phenotype for this category (for tissue area), 
               # then change Phenotype, Total Cells, Cell Density, others...
               # change all to 0 that are not consistent across all other rows of thiscss_tc
               new_css <- rbind(new_css, thiscss[alli_tc,]);
               lastrowi <- nrow(new_css);
               # Find all columns that are not the same with this Tissue Category, set to 0
               tf_notsame <- apply(thiscss_tc, 2, lengthUniqueNE1);
               new_css[lastrowi, tf_notsame] <- 0;
               new_css[lastrowi,'Phenotype'] <- ph;
               new_css[lastrowi,'Total Cells'] <- 0;
               new_css[lastrowi,'Cell Density (per megapixel)'] <- 0;
               new_css[lastrowi,'Nucleus Area (percent)'] <- "0.00%";
            }
         } # for each phenotype
      } # for each tissue category
   } # for each sample
   if(added_TC){
      # added a dummy tissue category as last column. Drop it before writing.
      new_css <- new_css[,1:(ncol(new_css)-1)];
   }
   #write.table(new_css, file="withZeros.txt", quote=FALSE, sep="\t", row.names=FALSE); 
   if(is.null(origfname)){
      outfile <- sprintf("%s_withZeros.txt", gsub(".txt$", "", basename(fname)));
   } else{
      outfile <- sprintf("%s_withZeros.txt", gsub(".txt$", "", origfname));
   }
   write.table(new_css, file=outfile, quote=FALSE, sep="\t", row.names=FALSE); 
   #oi <- order(new_css[,1], new_css[,3], new_css[,4], new_css[,5]);
   #write.table(new_css[oi,], file="withZeros_sorted.txt", quote=FALSE, sep="\t", row.names=FALSE); 
   #cell_seg_summary_per_patient(new_css);
   
   #return(new_css);
   return(sprintf("%s is ready.", outfile));
}



# Return true if all entries in vec are the same.
lengthUniqueEq1 <- function(vec){
   if(length(unique(vec)) == 1){
      return(TRUE);
   }
   return(FALSE);
}


# Return true if all entries in vec are the NOT the same.
lengthUniqueNE1 <- function(vec){
   if(length(unique(vec)) == 1){
      return(FALSE);
   }
   return(TRUE);
}


check_normal <- function(fname){
   #fname <- "DATA/MTC Panel 1 021618 Batch Merge_cell_seg_data.txt";
   csd <- read_cell_seg_data(fname); #  class(csd):   "tbl_df"     "tbl"        "data.frame"
   csd_mat <- as.matrix(csd);
   cn <- colnames(csd_mat);
   grepForNormal <- c("Normal", "NLN", "50/50", "NT", "normal");
   for (ci in c(1:5, 150)){
      unq <- unique(csd_mat[,ci]);
      print(sprintf("i=%i, u=%i  %s : %s", ci, length(unq), cn[ci], csd_mat[1,ci]), q=F);
      for(grep_string in grepForNormal){
         grepNormal <- grep(grep_string, unq, value=TRUE, fixed=TRUE);
         if(length(grepNormal) > 0){
            print(sprintf("grep_string: %s : n=%i : ", grep_string, length(grepNormal)), q=F);
            print(grepNormal);
         }
      }
   }
}

call_check_normal <- function(){
   #check_normal(fname = "DATA/120517 Panel 1 Batch Merge_cell_seg_data.txt");
   #check_normal(fname = "DATA/121217 Panel 1 Repeat Merge_cell_seg_data.txt");
   #check_normal(fname = "DATA/Panel 1 120707 merge_cell_seg_data.txt");
   check_normal(fname = "DATA/MTC Panel 1 021618 Batch Merge_cell_seg_data.txt");
   check_normal(fname = "DATA/MTC Panel 1 021618 Batch Merge_cell_seg_data.txt");
}



#options(shiny.maxRequestSize=30*1024^2); # 30MB for inputFile()
options(shiny.maxRequestSize=5000*1024^2); # 5GB for inputFile()

# Define UI ----
ui <- fluidPage(
  
  titlePanel("inForm2fcs"),

  sidebarLayout(
    sidebarPanel(h2("inForm2fcs"), 
    h4("Enable flow-style analysis of Vectra tissue scans for inForm segmentation results."),
    h4('Create fcs files, grouped by "Tissue Category" and "Phenotype" per each "Slide ID".')
# Read data produced by inForm software,
# Also want an "all" file per image. Including all phenotypes and tissue categories.
  ),
    #mainPanel("main panel",
    mainPanel("",
      #textOutput("tovalues"),
      #h1("First level title"),
      #h2("Second level title"),
      #h3("Third level title"),
      #h4("Fourth level title"),
      #h5("Fifth level title"),
      fileInput("file", h2("Select a file:"), accept=".txt", width="100%"),
      helpText(h3('Note: these files typically end in',
                    '"_cell_seg_data.txt"')),
# By default, Shiny limits file uploads to 5MB per file. You can modify this limit by using the shiny.maxRequestSize option. For example, adding options(shiny.maxRequestSize=30*1024^2) to the top of server.R would increase the limit to 30MB.
      #h6("Sixth level title"),
      br(),
      br(),
      br(),
      h3("Output will be sent to the current working directory:"),
      h3(code(getwd())),
      #h3('and will have the same file name as the input file'),
      #h3('with "_withZeros.txt" appended.'),
      br(),
      br(),
      h3('Please Note: This may take several minutes for large files (~9 minutes for a 1.5G file).'),
      h3('A summary will appear below when complete.'),
      #verbatimTextOutput("vtovalues")
      verbatimTextOutput("vtovalues", placeholder=TRUE)
      #textOutput("tovalues2")
      ) # mainPanel
  ) # sidebarLayout

  #You can use navbarPage to give your app a multi-page user interface that includes a navigation bar.
)


countRows <- function(fname){
   print("in countRows");
   #csd <- read_cell_seg_data(fname); #  class(csd):   "tbl_df"     "tbl"        "data.frame"
   csd <- read.table(fname, header=TRUE, sep="\t", as.is=TRUE, comment.char="", check.names=FALSE);
   return(nrow(csd));
}


# Define server logic ----
server <- function(input, output) {
   #output$tovalues <- renderPrint({
   #   req(input$file);
   #   str(input$file$datapath);
   #})
   output$vtovalues <- renderPrint({
   #output$vtovalues <- renderText({ 
      #req(input$file);
      #make_fcs(fname=input$file$datapath, origfname=input$file$name);
      if(is.null(input$file)){
         "A summary will appear here when complete."
      }else{
         make_fcs(fname=input$file$datapath, origfname=input$file$name);
      }
   }) 
}


# Run the app ---- from the command line.
shinyApp(ui = ui, server = server)
# start shiny app:  runApp("inFormA.R")


