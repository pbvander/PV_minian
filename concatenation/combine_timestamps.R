library(av)
library(tidyverse)
library(gtools)

######This code takes the "timeStamps.csv" file from multiple sessions that were concatenated together for analysis and returns one big timestamps file that will allow the output from MiniAn ran on concatenated sessions to be registered to "real" time and aligned to temperature data

##You should only need to run this code once, but you could also include a line in telemetry_miniscope analysis.R to source this code each time that that code is ran.

#Set working directory to [animal]/[date]/[session] (this folder should contain all the sessions that you want to concatenate)
setwd("C:/Users/paulv/Box/correalab/Member Folders/Paul Vander/Experiments/250417_circulating_E2_torpor_miniscope/pre-OVX_torpor/MT29/2025_05_22/session1")

##Miniscope timestamps
print("Miniscope timestamps")
ldf <- data.frame()
tsdf<- data.frame()
i=1
for (dir in list.dirs(recursive=F)[!grepl("concatenated", list.dirs(recursive=F))]%>%mixedsort()){
  print(dir)
  start_t = strsplit(dir, "/")[[1]][2]
  len=0
  for (avi in list.files(paste0(dir,"/My_V4_miniscope"))[endsWith(list.files(paste0(dir,"/My_V4_miniscope")),".avi")]%>%mixedsort()){
    len = len + av_media_info(paste0(dir,"/My_V4_miniscope/",avi))$video%>%pull(frames)
  }
  if (i==1){
    ldf[i, "start"] = 0
    ldf[i,"stop"] = len-1
    ldf[i, "start_time"] = start_t
    ldf[i, "length"] = len
  }
  if (i>1){
    ldf[i, "start"] = ldf%>%pull(length)%>%sum()
    ldf[i,"stop"] = ldf[i,"start"]+len-1
    ldf[i, "start_time"] = start_t
    ldf[i, "length"] = len
  }
  i=i+1
  ts <- read_csv(paste0(dir,"/My_V4_miniscope/","timeStamps.csv"),show_col_types = F)%>%
    rename(frame = `Frame Number`, time_ms = `Time Stamp (ms)`, buffer_index = `Buffer Index`)%>%mutate(start_time=start_t)
  tsdf<- rbind(tsdf, ts)
  if (len != ts%>%pull(frame)%>%max()+1){print("Warning! Timestamp max frame and length of frames in videos are not the same")}
}
tsdf<-merge(tsdf,ldf)%>%mutate(frame_og = frame, frame = frame_og + start)

write_csv(tsdf,"./concatenated/My_V4_miniscope/timeStamps.csv")
write_csv(ldf, "./concatenated/My_V4_miniscope/session_lengths.csv")

##Webcam timestamps
print("Webcam timestamps")
ldf <- data.frame()
tsdf<- data.frame()
i=1
for (dir in list.dirs(recursive=F)[!grepl("concatenated", list.dirs(recursive=F))]%>%mixedsort()){
  print(dir)
  start_t = strsplit(dir, "/")[[1]][2]
  len=0
  for (avi in list.files(paste0(dir,"/My_WebCam"))[endsWith(list.files(paste0(dir,"/My_WebCam")),".avi")]%>%mixedsort()){
    len = len + av_media_info(paste0(dir,"/My_WebCam/",avi))$video%>%pull(frames)
  }
  if (i==1){
    ldf[i, "start"] = 0
    ldf[i,"stop"] = len-1
    ldf[i, "start_time"] = start_t
    ldf[i, "length"] = len
  }
  if (i>1){
    ldf[i, "start"] = ldf%>%pull(length)%>%sum()
    ldf[i,"stop"] = ldf[i,"start"]+len-1
    ldf[i, "start_time"] = start_t
    ldf[i, "length"] = len
  }
  i=i+1
  ts <- read_csv(paste0(dir,"/My_WebCam/","timeStamps.csv"),show_col_types = F)%>%
    rename(frame = `Frame Number`, time_ms = `Time Stamp (ms)`, buffer_index = `Buffer Index`)%>%mutate(start_time=start_t)
  tsdf<- rbind(tsdf, ts)
  if (len != ts%>%pull(frame)%>%max()+1){print("Warning! Timestamp max frame and length of frames in videos are not the same")}
}
tsdf<-merge(tsdf,ldf)%>%mutate(frame_og = frame, frame = frame_og + start)

if (!dir.exists("./concatenated/My_WebCam")){dir.create("./concatenated/My_WebCam")} #Create My_WebCam directory, if it doesn't already exist

write_csv(tsdf,"./concatenated/My_WebCam/timeStamps.csv")
write_csv(ldf, "./concatenated/My_WebCam/session_lengths.csv")
