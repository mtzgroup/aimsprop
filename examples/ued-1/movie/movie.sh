ffmpeg -r 4 -i I%03d.png -c:v libx264 -pix_fmt yuv420p -r 10 out3.mp4
