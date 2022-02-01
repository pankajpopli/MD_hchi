start=1
end=200000
set xrange[-15:15]
set yrange[-15:15]
#set object rect from -8.5966,-7.4448 to 8.5966,7.4448
set si sq
set cbrange[0:3]
set cblabel "chi_i" font "helvetica, 18" 
do for [i=start:end:1]{
	set title  "Tri -> Rect, config:t =  ".i  font "helvatica,18"
	#set output "image.".i.".png"
	pl "time/new_cordinates_".i.".dat" u 2:3:8 w p pt 7 ps 2.0 palette notitle 
	#pause 0.0000001
   pause 0.05
}
