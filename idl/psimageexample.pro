pro psimageexample

image=bytscl(dist(512))			;create a dummy image

s=size(image)				;determine the pixel size
real_x = s(1)				;of the image
real_y = s(2)					

posbar=[0.1,.1,0.9,.3]			;set the position of the text region

posimage=[0.1,.4,0.9,0.9]			;set the position of the image

plot,position=posbar, $
	xstyle=4, $
	ystyle=4, $
	[0,1], $
	/nodata

devbar=fix(posbar*[!d.x_size, $
	!d.y_size, $
	!d.x_size, $
	!d.y_size])

xs=devbar(2)-devbar(0)-1
ys=devbar(3)-devbar(1)-1

plot,/noerase, $
	position=posbar, $
	xrange=[0,100], $
	yrange=[-1,1], $
	xstyle=1, $
	ystyle=1, $
	sin(indgen(100)/6.)

xyouts, .2,.4,"here is some sample text."
xyouts, .2,0.2,"it is drawn with xyouts in data coords."

devimage=fix(posimage*[!d.x_size, $
	!d.y_size, $
	!d.x_size, $
	!d.y_size])

xsi=devimage(2)-devimage(0)-1
ysi=devimage(3)-devimage(1)-1

plot,/noerase, $
	position=posimage, $
	xrange=[0,real_x], $
	yrange=[0,real_y], $
	xstyle=1, $
	ystyle=1, $
	/nodata, $
	[0,1]

if(!d.name eq 'PS') then begin
	tv,image,devimage(0)+1, $
	devimage(1)+1, $
	xsize=xsi, $
	ysize=ysi, $
	/device
endif else begin
	tv,congrid(image,xsi,ysi), $
	devimage(0)+1,devimage(1)+1
endelse
end


