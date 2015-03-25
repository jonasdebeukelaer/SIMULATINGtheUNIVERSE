cd ~/Documents/Uneyy/MSCI_PROJECT/;
python pmExecutable.py;
cd /Users/jonasdebeukelaer/Documents/Uneyy/MSCI_PROJECT/movies/SavingMovingMovie/;
mkdir movienGrid$1lBox$2; 
/Applications/VisIt.app/Contents/Resources/2.8.0/../bin/visit -cli -v 2.8.0 -reverse_launch -host 127.0.0.1 -port 5600 -key 475c4058c2c6776e443f ../../movieScript.py; 
ffmpeg -framerate 60 -i movie%04d.png -c:v libx264 -pix_fmt yuv420p movie.mp4;
open movie.mp4
