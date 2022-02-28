property idleCheck : 20 as integer
property idleCheck_usr : 120 as integer
set timer to 0
on idle
    --Check idle time
    set idletime to do shell script "ioreg -c IOHIDSystem | awk '/HIDIdleTime/ {print int($NF/1000000000); exit}'"
    set idletime to idletime as string
    set idletime to idletime as integer

    tell application "System Events"
        if idletime is less than idleCheck then (* 20 is 20 seconds. If a key was tapped within the idleCheck seconds, it quits the app. *)
            return idleCheck -- checks again in ... seconds
        else
            if idletime is greater than idleCheck_usr then (*  If a key was tapped after the idleCheck_usr seconds it opens the app. *)
                startCamera()
            end if

            return idleCheck
        end if
    end tell
end idle

-- run the iSenty program
on startCamera()
    tell application "Finder"
        activate

        open application file "iSentry.app" of folder "Applications" of startup disk
        delay 10
        tell application "System Events"
            tell process "iSentry"
                tell group 1 of window 1
                    tell button 1
                        click
                end tell
            end tell
        end tell
    end tell
end tell
end startCamera