#!/bin/sh

#  rsync_to_music_archive.sh
#  
#
#  Created by Jonathan Gordon on 3/26/15.

# check if shareddrive is connected
# use ControlPlane to handle IP connect http://www.controlplaneapp.com/

#handle mount here
#afp://Poopsmith._afpovertcp._tcp.local/Public/Shared Music
shareddrive="/Volumes/Public/Shared\ Music"


if [ ! -L $shareddrive ];
    then echo "=> drive not connected"
    echo "=> connecting ... making dir"
#osascript -e 'mount volume "smb://user:password@server/share"'
#mount_smbfs //user:password@server/share /tmp/mnt
echo "=> connected"
    else echo "=> drive connected"
fi

# check Music/archive folder
archivefolder=~/Music/archive/
cd $shareddrive
1 >> transfer.log




for f in $archivefolders/*; do
    album=$(basename $f)
    band=$(echo $album | awk -F "-" '{print $1}')
#if band folder exist copy to that folder if not create the folder
#find $shareddrive -iname "*$band*" -print
    if [ -d "$band" ];
        then rsync -azP "$f" "$band/"
        echo "Folder for $band exists: copied $album to $band `date +"%m-%d-%Y-%T"`" >> transfer.log
        else mkdir "$band";
        rsync -azP "$f" "$band/"
        echo "XX==> made folder: $band: copied $album to $band `date +"%m-%d-%Y-%T"`" >> transfer.log
    fi
done


shareddrive=~/Dropbox/test/;for f in $archivefolder/*; do album=$(basename $f); band=$(echo $album | awk -F "-" '{print $1}'); echo $band; query=$(find $shareddrive -iname "$band" -print ); echo $query; done