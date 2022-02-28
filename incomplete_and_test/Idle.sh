#!/bin/sh

# change for number of seconds to be defined as idle: now 10 min
VALUE=600
IDLE=$((`ioreg -c IOHIDSystem | sed -e '/HIDIdleTime/ !{ d' -e 't' -e '}' -e 's/.* = //g' -e 'q'` / 1000000000))


if [[ $IDLE -lt $VALUE ]]
    then
        #echo "not idle: active $IDLE seconds ago";
        exit 1;
    else
        #echo "Idle: inactive for $IDLE secs";
        exit 0;
fi