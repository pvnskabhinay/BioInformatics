library(eegkit)
eegcap()
eegcap(type = "2d")
data("eegcoord")
enames <- rownames(eegcoord)
length(enames) #87

eegcoord
write.csv(eegcoord,"C:\\Users\\pvnsk\\OneDrive\\Documents\\BI\\Project 2\\eegcoord.csv")
