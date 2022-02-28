#!/bin/sh

#  rsync_to_shared_drive.sh
#  
#
#  Created by Jonathan Gordon on 3/26/15.

# use ControlPlane to determine if connected to a work IP: http://www.controlplaneapp.com/
# and run this script automatically
# Actions > shellscript on arrival


# check if shareddrive is connected
shareddrive="/Volumes/steinlab"
if [ ! -L $shareddrive ];
    then echo "=> drive not connected"
    echo "=> connecting ... making dir"
    mkdir /Volumes/steinlab
    link=$COMIScred"@"$steindrive
    mount_smbfs //jargordo:D@ng3rs3ns3@med.uvm/shared/Labs/SteinLian/General/ /Volumes/steinlab
    echo "=> connected"
    else echo "=> drive already connected"
fi

# check papers folder (00_PAPER_REVISIONS) to see if there are links to Dropbox files
dbfolders='/Users/jgordon3/Dropbox/20*'
#only includes dated folder (20**_ prefix)

#for f in $dbfolders; do
#    subdir=$(basename $f)
#    papersfolder="/Users/jgordon3/WORK/00_PAPER_REVISIONS/$subdir"
#    if [ ! -L $papersfolder ]; then
#         echo "=> $subdir does not exist ... linking => $f to $papersfolder";
#        ln -s $f $papersfolder
#    fi
#done

#Sync terri paper
#sync log
#echo "-------------------------------------------------------" >> /Users/jgordon3/Dropbox/2017_Messier_drug_treatments/rsync_log.txt
#echo "=> synced steinlab <=> dropbox: " `date +"%m-%d-%Y-%T"` >> /Users/jgordon3/Dropbox/2017_Messier_drug_treatments/rsync_log.txt#
#echo "" >> /Users/jgordon3/Dropbox/2017_Messier_drug_treatments/rsync_log.txt

#rsync -rtuvv --copy-links /Users/jgordon3/Dropbox/2017_Messier_drug_treatments/* "/Volumes/steinlab/102 Manuscripts - Working Folder/0P.1. - Messier -Pfizer Drug Treatments" | grep -v 'uptodate' >> /Users/jgordon3/Dropbox/2017_Messier_drug_treatments/rsync_log.txt
#rsync -rtuvv --copy-links "/Volumes/steinlab/102 Manuscripts - Working Folder/0P.1. - Messier -Pfizer Drug Treatments/" /Users/jgordon3/Dropbox/2017_Messier_drug_treatments/ | grep -v 'uptodate' >> /Users/jgordon3/Dropbox/2017_Messier_drug_treatments/rsync_log.txt



#rsync WORK folder to shared drive
#rsync -azP --copy-links ~/WORK/* /Volumes/steinlab/Jonathan_Gordon/
echo "=> Last backup: " `date +"%m-%d-%Y-%T"` >> /Users/jgordon3/WORK/03_HPCC_VACC/logs/rsync_logs.txt
# sync to manuscipts
#append date to most recent files

#rsync -azuP --copy-links ~/Dropbox/2014_Wu_Histone_mods/* /Volumes/steinlab/102\ Manuscripts\ -\ Working\ Folder/Hai_Wu_BMSC_histone_2014/
#rsync -azuP --copy-links /Volumes/steinlab/102\ Manuscripts\ -\ Working\ Folder/Hai_Wu_BMSC_histone_2014/ ~/Dropbox/2014_Wu_Histone_mods/*
