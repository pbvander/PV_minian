##NEVER RUN THIS CODE ON ORIGINAL .AVI RECORDINGS!!!!!!!!!!!!!!!!!!!!!!
##This code is mean to combine miniscope videos from distinct directories into one for easy loading into MiniAn. This code irreversibly changes the name and location of any files (and will overwrite!!), so make sure to only run this code on a copy of the original .avi files. 
##It is recommended that you comment out the "file.rename" lines of this code at first, to test whether the assignment of new video names is working correct by inspecting printed output and concat_df.
##File structure: [animal]/[date]/[session]/concatenated/My_V4_Miniscope/[all .avi files numbered in order]

library(gtools)
library(dplyr)

#Set working directory (NOT PATH TO ORIGINAL VIDEOS - MAKE COPIES IN A NEW DIRECTORY)
setwd("C:/Users/paulv/Box/correalab/Member Folders/Paul Vander/Experiments/250417_circulating_E2_torpor_miniscope/pre-OVX_torpor/MT29/2025_05_22/session1/concatenated")

if ("concatenated" %in% strsplit(getwd(),"[/]")[[1]]){ #makes sure that code will only run if "conccatenated" is present in the current working directory path (to ensure that this code is only run on coapies of original files)

if (!dir.exists("./My_V4_Miniscope")){dir.create("./My_V4_Miniscope")}
  
#Combine files in one directory, renaming as needed
moved<-c()
i=1
ii=1
concat_df<-data.frame()
for (dir in list.dirs(recursive = F)%>%mixedsort()){
  for (file in list.files(paste0(dir,"/My_V4_Miniscope"))%>%mixedsort()){
    if (endsWith(file,".avi") & i==1){
      num<-strsplit(file,".",fixed=T)[[1]][1]
      print(paste0("Moving ",dir,"/My_V4_Miniscope/",file, " (no renaming)"))
      concat_df[ii,"original_name"]<-paste0(dir,"/",file)
      concat_df[ii,"new_name"]<-paste0("./",num,".avi")
      file.rename(from = paste0(dir,"/My_V4_Miniscope/",file),
                  to = paste0("./My_V4_Miniscope/",num,".avi"))
      moved<-c(moved,as.numeric(num))
      ii=ii+1
    }
    if (endsWith(file,".avi") & i>1){
      num<-strsplit(file,".",fixed=T)[[1]][1]
      print(paste0("Moving ",dir,"/My_V4_Miniscope/",file, " renaming to ",max(moved)+1,".avi"))
      concat_df[ii,"original_name"]<-paste0(dir,"/",file)
      concat_df[ii,"new_name"]<-paste0("./",max(moved)+1,".avi")
      file.rename(from = paste0(dir,"/My_V4_Miniscope/",file),
                  to = paste0("./My_V4_Miniscope/",max(moved)+1,".avi"))
      moved<-c(moved,max(moved)+1)
      ii=ii+1
    }
  }
  i=i+1
}

write.csv(concat_df,"./concatenation code output.csv", row.names = F)
}
