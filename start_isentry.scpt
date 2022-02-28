-- on startCamera()
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
-- end startCamera