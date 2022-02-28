#!/bin/sh

# monitor_folder.sh
# script monitors a path ($path) and looks for new data and message to Slack
# requires inotify-tools package: sudo apt-get install inotify-tools
# Created by Jonathan Gordon on 6/27/18.

#monitored path
path="/slipstream/home/agtc"
#path="/slipstream/home/jonathan"


inotifywait -m $path -e create -e moved_to |
 while read path action file; do
 message="New data: ${file} appeared in the directory ${path}. Do something with it or don't... I don't care... I'm a bot"
 eval $(echo 'curl -X POST -H 'Content-type: application/json' --data "{\"text\":\"${message}\"}" https://hooks.slack.com/services/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')
done

