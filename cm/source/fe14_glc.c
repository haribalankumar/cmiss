
#include <gl/gl.h>
#include <gl/device.h>
#include <math.h>

/* Global declarations */

#define CURRENT_TRANS 0
#define MAX_PICKS 20
#define MAX_LOCATORS 20
#define MAX_CHOICES 20
#define MAX_POLYLINES 256
#define MAX_POLYMARKERS 256
#define MAX_WINDOWS 35
#define MAX_SEGMENTS 500
#define MAX_TEXTS 256
#define MAX_TIMES 50  /* for animation */

int multi_mode=0;
int qdevices=0;
short mouse_position[2];
Device mouse[2]={MOUSEX,MOUSEY};
short iseg[MAX_SEGMENTS]={0}; /*keeps track of current graphical objects,
	           their detectablility and visibility
		   iseg[iw][isegnum]=0 for not defined, 1 for invisible, 2 for
		   visible and 3 for detectable*/
short iwseg[MAX_SEGMENTS]={0}; /*keeps track of windows for objects*/
short iseg_time[MAX_TIMES][MAX_SEGMENTS]={0}; /*for animation */
short max_segments[MAX_TIMES]={0};            /*for animation*/
int frame=0;            /*animation frame number*/
int animate=0;          /*animation flag*/

struct Gl{
   int idws[MAX_WINDOWS];             /*gl identifiers of open windows*/
   int background_color[MAX_WINDOWS]; /*colour indices of window backgrounds*/
   char *titlews[MAX_WINDOWS];        /*titles of each window*/
   Matrix modelview[MAX_WINDOWS];     /*model/viewing transformation*/
   float wc_window[MAX_WINDOWS][7];   /*world coord window for each ws */
   short view_port[MAX_WINDOWS][4];     /*viewport in pixel coordinates*/
   float diagonal;                    /*length of diagonal of wc window*/
   float centre[MAX_WINDOWS][3];    /*viewing point*/
   short colourmap[MAX_WINDOWS][256][3]; /*colour luts */
}gl={0};

struct Polyline_rep{
   short type[MAX_POLYLINES];
   float width[MAX_POLYLINES];
   unsigned short colour[MAX_POLYLINES];
}polyline_rep[MAX_WINDOWS];

struct Polymarker_rep{
   short type[MAX_POLYMARKERS];
   float size[MAX_POLYMARKERS];
   unsigned short colour[MAX_POLYMARKERS];
}polymarker_rep[MAX_WINDOWS];

struct Text_rep{
   short font[MAX_TEXTS];
   float width[MAX_TEXTS];
   float spacing[MAX_TEXTS];
   unsigned short colour[MAX_TEXTS];
}text_rep[MAX_WINDOWS];

struct Pick{
   int number_active;
   int window_id[MAX_PICKS];
}pick_devices={0};

struct Locator{
   int number_active;
   int window_id[MAX_LOCATORS];
}locator_devices = {0};

struct Choice{
   int number_active;
   int window_id[MAX_CHOICES];
}choice_devices = {0};

unsigned short font_bits[] = {
   /* + */
   0x2000,0x2000,0xf800,0x2000,0x2000
};

#define font_nr (sizeof font_bits)
#define ASSIGN(fontch, of, wi, he, xof, yof, wid) \
	       fontch.offset = of; \
	       fontch.w = wi; \
	       fontch.h = he; \
	       fontch.xoff = xof; \
	       fontch.yoff = yof; \
	       fontch.width = wid

Fontchar font_chars[127];

foreground_()
{
  foreground();
}


fbclear_(iw)

   int *iw;

{
  color(0);
  clear();
  color(1);
}


fbpl_(iw,npt,xpl,ypl)

   int *iw,*npt;
   float *xpl,*ypl;
{
   float points[200][3];
   int i;

   for (i=0; i< *npt && i<200; i++ ) {
     points[i][0] = *(xpl+i);
     points[i][1] = *(ypl+i);
     points[i][2] = 0;
   }
   gl_polyline_(npt,points);
}


fbpm_(iw,npt,xpl,ypl)
 
   int *iw,*npt;
   float *xpl,*ypl;
{
   float points[200][3];
   int i;
 
   for (i=0; i< *npt && i<200; i++ ) {
     points[i][0] = *(xpl+i);
     points[i][1] = *(ypl+i);
     points[i][2] = 0;
   }  
   gl_polymarker_(npt,points);
}
     

gl_acwk_(iw)          
  
   /* activate window */
   int *iw; /* workstation window number iw=1,2,3*/

{
   winset( gl.idws[*iw] );
   /*setmap(*iw);*/
   /* set appropriate transformation matrix here? */
   gl_select_xform_(iw);

   makeobj(CURRENT_TRANS);
   viewport(gl.view_port[*iw][0],gl.view_port[*iw][1],
   gl.view_port[*iw][2],gl.view_port[*iw][3]);
   gl_select_xform_(iw);
   closeobj();

}


gl_animate_(iw,ntimes,times,delay,loop)

   /* use mouse to update viewing transformation */
   int *iw,*ntimes,*times,*delay,*loop;

{
   int cont,dev,x,y,time;
   short val;
   float factor;

   frontbuffer(FALSE);
   swapinterval((short) *delay);
   /* event loop */
   cont = 1;
   time = 0;
   while (cont) {
     /* draw at this time */
     gl_draw_struct_at_time_(iw,times+time);
     swapbuffers();
     if (qtest()) {
       dev = qread(&val);
       if (dev==LEFTMOUSE) {
         swapinterval(0);
	 /*get current cursor location*/
	 x = getvaluator(MOUSEX);

	 /*if z key is pressed zoom*/
	 if( getbutton(ZKEY) ) {
           /*while left button down*/
           while (getbutton(LEFTMOUSE)) {
	     /*poll mouse position and zoom*/
             gl_select_xform_(iw);
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     factor = getvaluator(MOUSEX)-x;
	     if( factor >= 0.0 )
	       factor = 0.01*factor+1;
	     else
	       factor = 100.0/(100.0 - factor);
	     scale(factor,factor,factor);
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
             gl_draw_struct_at_time_(iw,times+time);
	     swapbuffers();
           } 
	 }
	 /*if t key is pressed translate*/
         else if( getbutton(TKEY) ) {
	   factor = 0.0;
           /*while left button down*/
           while (getbutton(LEFTMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
	     factor = (getvaluator(MOUSEX)-x)*gl.diagonal*0.001;
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     translate( factor,0.0,0.0);
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
             gl_draw_struct_at_time_(iw,times+time);
	     swapbuffers();
           } 
           gl.centre[*iw][0] = gl.centre[*iw][0] - factor;
	 }
	 else { /*rotate*/
           /*while left button down*/
           while (getbutton(LEFTMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     rot(0.3*(getvaluator(MOUSEX)-x), 'x');
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
             gl_draw_struct_at_time_(iw,times+time);
	     swapbuffers();
           } 
	 }
	 getmatrix(gl.modelview[*iw] );
         swapinterval((short) *delay);
       }
       else if (dev==MIDDLEMOUSE) {
         swapinterval(0);
	 x = getvaluator(MOUSEX);
	 /*if t key is pressed translate*/
         if( getbutton(TKEY) ) {
	   factor = 0.0;
           /*while left button down*/
           while (getbutton(MIDDLEMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
	     factor = (getvaluator(MOUSEX)-x)*gl.diagonal*0.001;
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     translate( 0.0,factor,0.0);
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
             gl_draw_struct_at_time_(iw,times+time);
	     swapbuffers();
           } 
           gl.centre[*iw][1] = gl.centre[*iw][1] - factor;
	 }
	 else { /*rotate*/
           /*while left button down*/
	   while (getbutton(MIDDLEMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     rot(0.3*(getvaluator(MOUSEX)-x), 'y');
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
             gl_draw_struct_at_time_(iw,times+time);
	     swapbuffers();
           } 
	 }
	 getmatrix(gl.modelview[*iw] );
         swapinterval((short) *delay);
       }
       else if (dev==RIGHTMOUSE) {
         swapinterval(0);
	 x = getvaluator(MOUSEX);
	 /*if t key is pressed translate*/
         if( getbutton(TKEY) ) {
	   factor = 0.0;
           /*while left button down*/
           while (getbutton(RIGHTMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
	     factor = (getvaluator(MOUSEX)-x)*gl.diagonal*0.001;
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     translate( 0.0,0.0,factor);
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
             gl_draw_struct_at_time_(iw,times+time);
	     swapbuffers();
           } 
           gl.centre[*iw][2] = gl.centre[*iw][2] - factor;
	 }
	 else { /*rotate*/
           /*while left button down*/
	   while (getbutton(RIGHTMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     rot(0.3*(getvaluator(MOUSEX)-x), 'z');
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
             gl_draw_struct_at_time_(iw,times+time);
	     swapbuffers();
           } 
	 }
	 getmatrix(gl.modelview[*iw] );
         swapinterval((short) *delay);
       }
       else if (dev==NKEY && getbutton(NKEY)) {
	 time=time+1;
	 if(time == *ntimes)time = 0;
       }
       else if (dev==ESCKEY && getbutton(ESCKEY)) cont = 0;
     }   
     if(*loop==1)time = time+1;
     if(time == *ntimes)time = 0;
   } 
   frontbuffer(TRUE);
   swapinterval(0);
   depthcue(FALSE);
}


gl_build_xform_matrix3_(fixed_pt,shift,angle1,angle2,angle3,scale,status,matrix)

  /* put shift, rotate and scale into one 4x4 matrix */
  float *fixed_pt;               /* not used */
  float *shift;                  /* array of 3 shifts */
  float *angle1,*angle2,*angle3; /* 3 rotation angles */
  float *scale;                  /* array of 3 shifts */
  int *status;                   /* error flag */
  float *matrix;                 /* 4x4 transform matrix */

{

}


gl_cell_array_(xmin_ca,ymin_ca,xmax_ca,ymax_ca,nxm,nym,nxdim,nydim,image)
 
   /* draw image in i2p */
   int *xmin_ca,*ymin_ca,*xmax_ca,*ymax_ca,*nxm,*nym,*nxdim,*nydim,*image;

{
   int i,shift8,shift16,mask;

   cmode();
   drawmode(NORMALDRAW);
   gconfig();
   mask = 255;
   /*for( i=0;i<512*512; i++) {
     shift8 = *(image+i)  << 8;
     shift16 = *(image+i)  << 16;
     *(image+i) = *(image+i) + shift8 + shift16;
   }*/
   for( i=0;i<512*512; i++) *(image+i) = *(image+i) + 512;
   lrectwrite(0,0,511,511,(unsigned long *)image); 
   /*for( i=0;i<512*512; i++) *(image+i) = *(image+i) & mask;*/
   for( i=0;i<512*512; i++) *(image+i) = *(image+i) - 512;
   drawmode(OVERDRAW);
   gconfig();
   mapcolor(1,255,0,0);
   color(1);
}


gl_clear_ws_(iw)

   /* clear window */
   int *iw;  /* window identifier */

{
   int index;

   winset(gl.idws[*iw]);
   index = (gl.background_color[*iw]);
   c3s( gl.colourmap[*iw][index] );
   clear();
}


gl_close_segment_(isegnum)

   /* close segment */
   int *isegnum;  /* segment identifier */

{
   closeobj();
   callobj(*isegnum);
}


gl_close_ws_(iw)

   /* close window */
   int *iw;  /* window identifier */

{
   winclose(gl.idws[*iw]);
}


gl_compose_matrix3_(m1,m2,status,mout)

   /*premultiply m2 with m1 result to mout*/
   float *m1,*m2,*mout;
   int *status;

{
   /*may not need this if it is just called from set view */
}


gl_create_segment_(isegnum,iw)

   /* create segment */
   int *isegnum;  /* segment identifier */
   int *iw;

{
   makeobj(*isegnum);
   loadname((short)*isegnum);
   if(*isegnum>MAX_SEGMENTS-1)
     printf(" too many segments");
   else {
     iseg[*isegnum]=2;
     iwseg[*isegnum] = *iw;
   }
   if(animate) {
     iseg_time[frame][max_segments[frame]] = *isegnum;
     max_segments[frame]++;
   }
   
}


gl_dawk_(iw)          
  
   /* deactivate window */
   int *iw; /* workstation window number iw=1,2,3*/

{
   /*seems gl can only have one window active at a time*/
   /*so this is not supported*/
}


gl_delete_choice_(iw)

   /* delete choice menu */
   int   *iw; /* workstation window */

{
   int i,j;

   /* search for this window id*/
   for (i=0;i<choice_devices.number_active;i++) {
     if (choice_devices.window_id[i]==*iw) {
       /*pull rest down*/
       for (j=i;j<choice_devices.number_active-1;j++) 
	 choice_devices.window_id[j]=choice_devices.window_id[j+1]; 
       choice_devices.number_active--;
     }
   }

}


gl_delete_locator_(iw)

   /* delete locator */
   int   *iw; /* workstation window */

{
   int i,j;

   /* search for this window id*/
   for (i=0;i<locator_devices.number_active;i++) {
     if (locator_devices.window_id[i]==*iw) {
       /*pull rest down*/
       for (j=i;j<locator_devices.number_active-1;j++) 
	 locator_devices.window_id[j]=locator_devices.window_id[j+1]; 
       locator_devices.number_active--;
     }
   }
}


gl_delete_pick_(iw)

   /* delete pick */
   int   *iw; /* workstation window */

{
   int i,j;

   /* search for this window id*/
   for (i=0;i<pick_devices.number_active;i++) {
     if (pick_devices.window_id[i]==*iw) {
       /*pull rest down*/
       for (j=i;j<pick_devices.number_active-1;j++) 
	 pick_devices.window_id[j]=pick_devices.window_id[j+1]; 
       pick_devices.number_active--;
     }
   }
}


gl_delete_segment_(isegnum,iw)

   /* delete segment */
   int *isegnum; /* segment identifier */
   int *iw;

{
   delobj(*isegnum);
   iseg[*isegnum]=0;
}


gl_delete_stroke_(iw)

   /* delete stroke */
   int   *iw; /* workstation window */

{

}


gl_delete_valuator_(iw)

   /* delete valuator window */
   int   *iw; /* workstation window */

{

}


gl_detect_(isegnum,detect,iw)

   /* delete valuator window */
   int   *isegnum,*detect,*iw; 

{
   iseg[*isegnum] = *detect;
}


gl_eval_view_orien_matrix3_(view_ref_pt,view_plane,view_up,status,matrix)

   /* calculate the viewing transformation */
   float *view_ref_pt; /*is the point in the object at which you look*/ 
   float *view_plane;  /*defines the normal to the plane at which you look*/
   float *view_up;     /*which way is up*/
   int *status;        /*error flag*/
   float *matrix;      /*viewing matrix (4x4)*/

{

}


gl_eval_view_map_matrix3_(window,viewport,mode_proj,proj_ref_pt,
   view_plane_dist,back_plane_dist,front_plane_dist,status,matrix)

   /* calculate the mapping projection */
   float *window; /*area of the view plane which gets mapped onto the viewport*/ 
   float *viewport;/*area of workstation window in which the object will appear*/
   int *mode_proj; /*orthographic or perspective*/
   float *proj_ref_pt; /*the point from which you look*/
   float *view_plane_dist,*back_plane_dist,*front_plane_dist;
   int *status;    /*error flag*/
   float *matrix;  /*mapping matrix (4x4)*/

{
   /* may need to merge this with gl_eval_view_orien_matrix3 */
}


gl_event_(iw,class,idevice,status)

   /*get next event from event queue if present */
   int  *iw;      /* window identifier for input device */
   int  *class;   /* class of input device:
		     0:none
		     1:choice
		     2:locator
		     3:pick
		     4:string
		     5:valuator*/
   int  *idevice; /* type of input device */
   int  *status;  /* 0 for middlebutton (end input condition), 1 otherwise*/

{
   short looking,device,val;

   /*test event queue*/
   *status=1;
   *class=0;
   looking=1;
   while (looking) {      /*keeping looking until queue empty or button down*/
     device=qtest();
     /* return immediately if no event present*/
     if(device==0) return;
     if(device == LEFTMOUSE || device == RIGHTMOUSE || device == MIDDLEMOUSE){
       device=qread(&val);/*get status of button*/
       if(val)            /*if button is down stop looking otherwise continue*/
         looking=0;
     }
     else
       qread(&val);
   }

   switch (device) {
   case LEFTMOUSE:
     if (val==1) {
       getdev( 2, mouse, mouse_position );
       /*check against open devices*/
       if ( in_choice(mouse_position,iw) ) 
	 *class=1;
       else if ( in_locator(mouse_position,iw) ) 
	 *class=2;
       else if ( in_pick(mouse_position,iw) )
	 *class=3;
       else
	 *class=0;
     }
     break;
   case MIDDLEMOUSE:
     if(val==1) {
       getdev( 2, mouse, mouse_position );
       if ( in_choice(mouse_position,iw) ) 
	 *class=1;
       else
         /*status=0 is flag for end of input*/
         *status=0;
     }
     break;
   case RIGHTMOUSE:
     if (val==1) {
       getdev( 2, mouse, mouse_position );
       /*check against open devices*/
       if ( in_choice(mouse_position,iw) ) 
	 *class=1;
       else if ( in_locator(mouse_position,iw) ) 
	 *class=2;
       else if ( in_pick(mouse_position,iw) )
	 *class=3;
       else
	 *class=0;
     }
     break;
   }

}


gl_fill_area_(npoints,points)

   /* draw filled area */
   int npoints;    /* number of points in points */
   float *points;  /* array of points */

{

}


gl_flush_device_events_(class,flag)

  /* reset event queue */
  char *class;  /* input class */
  int *flag;    /* ? */

{

}


gl_get_choice_(iw,noch,status)

   /*get choice option number */
   int *iw;     /* window identifier */
   int *noch;   /* option number */
   int *status; /* 1 for input ok, 0 for error */

{
   /* an event may have happened which was a choice in which case the x,y */
   /* location of the choice needs to be checked against option areas */
   *noch = 0;
   *status = 0;
   if(getbutton(LEFTMOUSE)) 
     *noch = 1;
   if(getbutton(MIDDLEMOUSE)) 
     *noch = 2;
   if(getbutton(RIGHTMOUSE)) 
     *noch = 3;
   if(*noch>0) *status = 1;
 
}


gl_get_locator_(iw,x,y,status)

   /*get locator position */
   int *iw;       /* window identifier */
   float *x,*y;   /* location in world coords? */
   int *status; /* 1 for input ok, 0 for error */

{
   float wx1,wy1,wz1,wx2,wy2,wz2;

   /* an event may have happened which was a locator in which case the x,y */
   /* location needs to be checked against iw area */
   getdev( 2, mouse, mouse_position );
   mapw(CURRENT_TRANS,mouse_position[0],mouse_position[1],
     &wx1,&wy1,&wz1,&wx2,&wy2,&wz2);
   *x=wx1;
   *y=wy1;
}


gl_get_pick_(iw,id,status,object,pickid)

   /*get picked object */
   int *iw;     /* window identifier */
   int *object; /* object number */
   int *id;     /* object id */
   int *pickid; /* pick id */
   int *status; /* 1 for input ok, 0 for error */

{
   /* an event may have happened which was a pick in which case the x,y*/
   /* location of the cursor needs to be checked against object areas */
   int numpicks;
   short buffer[20],i,j,items,name,index;

   /*zero the name stack*/
   initnames();
   /*get into picking mode*/
   pick(buffer,20);
   /*callobj(CURRENT_TRANS);*/
   gl_select_xform_(iw);

   /*draw all segments*/
   for (i=1;i<MAX_SEGMENTS-1;i++) {
     if( isobj(i) ) 
       callobj(i);
   }
   /*get out of picking mode*/
   numpicks=endpick(buffer);
   
   /*loop over picked namelists*/
   index=0;
   *object=0; 
   for( i=0;i<numpicks;i++) {
     items=buffer[index++];
     for (j=0;j<items;j++) {
       name=buffer[index++];
       if( iseg[name]==3 /*detectable*/) {
	 *object=name;
	 return;
       }
     }
   }
}


gl_get_string_(iw,string,size,status)

   /*get string */
   int *iw;      /* window identifier */
   char *string; /* string */
   int *size;    /* number of characters */
   int *status; /* 1 for input ok, 0 for error */

{
   /* dinner isn't ready yet*/
}


gl_get_valuator_(iw,value,status)

   /*get value */
   int *iw;      /* window identifier */
   float *value; /* value */
   int *status; /* 1 for input ok, 0 for error */

{
   /* an event may have happened which was a valuator in which case the x,y*/
   /* location of the cursor needs to be checked against valuator areas */
   /* In sample mode just get the current cursor location */
}


gl_init_choice_(iw,idevice,ecarea)

   /* initialise choice menu */
   int   *iw;       /* workstation window */
   int   *idevice;  /* device type */
   float *ecarea;   /* screen coord window area */ 

{
   choice_devices.number_active++;
   if (choice_devices.number_active>MAX_LOCATORS) {
     printf(" too many active choices");
     choice_devices.number_active=MAX_CHOICES;
     return;
   }
   choice_devices.window_id[choice_devices.number_active-1] = *iw;

}

gl_init_locator_(iw,idevice,xref,yref,x1,x2,y1,y2)

   /* initialise locator */
   int   *iw;               /* workstation window */
   int   *idevice;          /* device type */
   float *xref,*yref;       /* screen coord initial point */ 
   float *x1,*y1,*x2,*y2;   /* screen coord window limits */ 

{
   locator_devices.number_active++;
   if (locator_devices.number_active>MAX_LOCATORS) {
     printf(" too many active locators");
     locator_devices.number_active=MAX_LOCATORS;
     return;
   }
   locator_devices.window_id[locator_devices.number_active-1] = *iw;

}


gl_init_pick_(iw,idevice)

   /* initialise pick */
   int   *iw;       /* workstation window */
   int   *idevice;  /* device type */

{
   pick_devices.number_active++;
   if (pick_devices.number_active>MAX_PICKS) {
     printf(" too many active picks");
     pick_devices.number_active=MAX_PICKS;
     return;
   }
   pick_devices.window_id[pick_devices.number_active-1] = *iw;
}


gl_init_stroke_(iw,idevice)

   /* initialise stroke */
   int   *iw;       /* workstation window */
   int   *idevice;  /* device type */

{

}


gl_init_valuator_(iw,idevice,ecarea,valmin,valmax)

   /* initialise valuator */
   int   *iw;       /* workstation window */
   int   *idevice;  /* device type */
   float *ecarea;   /* screen coord window area */ 
   float *valmin,*valmax;  /* limits of values returned */ 

{

}


gl_init_view_(iw)

   /*initialise viewing transformations, world coord window*/
   int *iw;

{
   float vrpx,vrpy,vrpz,prpx,prpy,prpz,xmin,xmax,ymin,ymax,zmin,zmax;
   
   vrpx = 0.5*(gl.wc_window[*iw][0]+gl.wc_window[*iw][1]);
   vrpy = 0.5*(gl.wc_window[*iw][2]+gl.wc_window[*iw][3]);
   vrpz = 0.5*(gl.wc_window[*iw][4]+gl.wc_window[*iw][5]);
   prpx = vrpx;
   prpy = vrpy;
   prpz = 3*gl.wc_window[*iw][6] + vrpz;
   xmin = gl.wc_window[*iw][0] - vrpx;
   xmax = gl.wc_window[*iw][1] - vrpx;
   ymin = gl.wc_window[*iw][2] - vrpy;
   ymax = gl.wc_window[*iw][3] - vrpy;
   zmin = gl.wc_window[*iw][6] * 2.5;
   zmax = gl.wc_window[*iw][6] * 3.5;
   ortho( xmin,xmax,ymin,ymax,zmin,zmax);
   lookat( prpx,prpy,prpz,vrpx,vrpy,vrpz,0.0);
   getmatrix( gl.modelview[*iw] );
}


in_locator(mouse_position,iw)

   /*check mouse position with active locators*/
   short mouse_position[];
   int *iw;

{
   int i,id_ws;

   /*loop over active locators*/
   for (i=locator_devices.number_active;i>0;i--) {
     /*check if in viewport*/
     id_ws=locator_devices.window_id[i-1];
     if (gl.view_port[id_ws][0]<=mouse_position[0] &&
	 gl.view_port[id_ws][1]>=mouse_position[0] &&
	 gl.view_port[id_ws][2]<=mouse_position[1] &&
	 gl.view_port[id_ws][3]>=mouse_position[1] ) {
       *iw=id_ws;
       return 1;
     }
   }
   return 0;
}


in_pick(mouse_position,iw)

   /*check mouse position with active locators*/
   short mouse_position[];
   int *iw;

{
   int i,id_ws;

   /*loop over active pick*/
   for (i=pick_devices.number_active;i>0;i--) {
     /*check if in viewport*/
     id_ws=pick_devices.window_id[i-1];
     if (gl.view_port[id_ws][0]<=mouse_position[0] &&
	 gl.view_port[id_ws][1]>=mouse_position[0] &&
	 gl.view_port[id_ws][2]<=mouse_position[1] &&
	 gl.view_port[id_ws][3]>=mouse_position[1] ) {
       *iw=id_ws;
       return 1;
     }
   }
   return 0;
}


in_choice(mouse_position,iw)

   /*check mouse position with active choice devices*/
   short mouse_position[];
   int *iw;

{
   int i,id_ws;

   /*loop over active choice*/
   for (i=choice_devices.number_active;i>0;i--) {
     /*check if in viewport*/
     id_ws=choice_devices.window_id[i-1];
     if (gl.view_port[id_ws][0]<=mouse_position[0] &&
	 gl.view_port[id_ws][1]>=mouse_position[0] &&
	 gl.view_port[id_ws][2]<=mouse_position[1] &&
	 gl.view_port[id_ws][3]>=mouse_position[1] ) {
       *iw=id_ws;
       return 1;
     }
   }
   return 0;
}



gl_inq_colour_(status,ncolours,flag,nindices)

   /* inquire number of colour indices in LUT */
   int *status,*ncolours,*flag,*nindices;

{
   /* GBS changed from GD_BITS_NORM_SNG_MMAP for IBM compat */
   *nindices=pow(2,getgdesc(GD_BITS_NORM_SNG_CMODE));
   *status=0;
}


gl_inq_colour_rep_(iw,colour_index,status,col)

   /* inquire rgb representation of index entry in LUT */
   int *iw,*colour_index,*status;
   float *col;  

{
   Colorindex index;
   short r,g,b;
   /* set the current colour lut here */
   r = gl.colourmap[*iw][*colour_index][0];
   g = gl.colourmap[*iw][*colour_index][1];
   b = gl.colourmap[*iw][*colour_index][2];
   *col=(float)r;
   *(col++)=(float)g;
   *(col++)=(float)b;
}


gl_inq_disp_(status,xdisp,ydisp)

   /* inquire size of display screen */
   int *status;
   float *xdisp,*ydisp; /* size of screen */

{
   *xdisp=(float)getgdesc(GD_XPMAX);
   *ydisp=(float)getgdesc(GD_YPMAX);
   *status=1;
}


gl_inq_elem_(root,elem,status,type,size)

   /* inquire nature of graphical object */
   int *root,*elem,*status,*type,*size;

{

}


gl_inq_ws_(iw,status)

   /* inquire status of iw */
   int *iw,*status;

{

}


gl_invis_(iw,index,status)

   /* set invis status of index*/
   int *iw,*index,*status;

{
   iseg[*index] = *status ;
   gl_redraw_all_struct_(iw);
}


gl_open_(ioerr)

   /* open graphics */
   int *ioerr; /* file unit for error messages */

{
   int i;

   /* initialise titles */
   for (i=0;i<MAX_WINDOWS;i++) {
     gl.titlews[i]=" cmiss 1";
   }
   pick_devices.number_active=0;
   locator_devices.number_active=0;
}


gl_open_wkst_(iw)

   /* open workstation window */
   int *iw; /* window identifier */

{
   int index;
   static short grey[]={0,0,100};

   gl.idws[*iw]=winopen(gl.titlews[*iw]);
   mmode(MVIEWING);
   /*multimap();*/
   zbuffer(TRUE); 
   zclear();
   RGBmode(); /*is not compatible with multimap!*/
   doublebuffer();  /* might need to put this with multimap*/
   gconfig();
   frontbuffer(TRUE);
   multi_mode=1;
   for(index=0;index<256;index++) 
     mapcolor(index+512,index,index,index);
   if( qdevices==0 ) {
     qdevice(LEFTMOUSE);
     qdevice(MIDDLEMOUSE);
     qdevice(RIGHTMOUSE);
     qdevice(ESCKEY);
     qdevice(NKEY);
     ASSIGN(font_chars['+'], 0, 5, 5, -2, -2, 5);
     defrasterfont(1,5,127,font_chars,font_nr,font_bits);
     qdevices=1;
   }
   gl.background_color[*iw] = 0;
   index = (gl.background_color[*iw]);
   c3s( gl.colourmap[*iw][index] );
   clear();
}


gl_polyline_(npoints,points)

   /*draw line connecting points */
   int *npoints;  /* number of points in points */
   float *points; /* array of points */

{
   int i;
   /*assume points are stored with x,y,z in sequential order*/
   bgnline();
      for (i=0;i<*npoints;i++) 
	 v3f(points+(i*3));
   endline();
}


gl_polymarker_(npoints,points)

   /*draw markers at points */
   int *npoints;   /* number of points in points */
   float *points;  /* array of points */

{
   int i;
   float p[3],q[3];

   /*assume points are stored with x,y,z in sequential order*/
   font(1);
   for (i=0;i<*npoints;i++) {
     cmov(*(points+(i*3)),*(points+(i*3)+1),*(points+(i*3)+2));
     charstr("+");
   }
   font(0);

     /*draw three lines each 0.01*size of window
     p[0]= *(points+(i*3)) - 0.005*(gl.diagonal);
     p[1]= *(points+(i*3)+1);
     p[2]= *(points+(i*3)+2);
     q[0]= *(points+(i*3)) + 0.005*(gl.diagonal);
     q[1]= *(points+(i*3)+1);
     q[2]= *(points+(i*3)+2);
     bgnline();
       v3f(p);
       v3f(q);
     endline();
     p[0]= *(points+(i*3));
     p[1]= *(points+(i*3)+1) - 0.005*(gl.diagonal);
     p[2]= *(points+(i*3)+2);
     q[0]= *(points+(i*3));
     q[1]= *(points+(i*3)+1) + 0.005*(gl.diagonal);
     q[2]= *(points+(i*3)+2);
     bgnline();
       v3f(p);
       v3f(q);
     endline();
   }*/

}


gl_draw_struct_at_time_(iw,time)

   /*redraw everything at time for animations */
   int *iw,*time;

{
   int i;
   int index;

   /* clear viewport and redraw all segments*/
   index = (gl.background_color[*iw]);
   c3s( gl.colourmap[*iw][index] );
   clear();
   zclear();
   for (i=0;i<max_segments[*time];i++) {
     if( isobj(iseg_time[*time][i]) ) 
       callobj(iseg_time[*time][i]);
   }
}





gl_redraw_all_struct_(iw)

   /* redraw everything on iw */
   int *iw;

{
   int i;
   int index;

   /* clear viewport and redraw all segments*/
   index = (gl.background_color[*iw]);
   c3s( gl.colourmap[*iw][index] );
   clear();
   zclear();
   for (i=1;i<MAX_SEGMENTS-1;i++) {
     if( isobj(i) && iseg[i]==2 && iwseg[i] == *iw ) 
       callobj(i);
   }
}


gl_rotate_(iw)

   /* use mouse to update viewing transformation */
   int *iw;

{
   int cont,dev,x,y;
   short val;
   float factor;

   frontbuffer(FALSE);
   /*
   lsetdepth(0,0x7f0000);
   lRGBrange(0,0,0,255,255,255,0,0x7f0000);
   depthcue(TRUE);
   */
   /* event loop */
   cont = 1;
   while (cont) {
     if (qtest()) {
       dev = qread(&val);
       if (dev==LEFTMOUSE) {
	 /*get current cursor location*/
	 x = getvaluator(MOUSEX);

	 /*if z key is pressed zoom*/
	 if( getbutton(ZKEY) ) {
           /*while left button down*/
           while (getbutton(LEFTMOUSE)) {
	     /*poll mouse position and zoom*/
             gl_select_xform_(iw);
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     factor = getvaluator(MOUSEX)-x;
	     if( factor >= 0.0 )
	       factor = 0.01*factor+1;
	     else
	       factor = 100.0/(100.0 - factor);
	     scale(factor,factor,factor);
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
  	     gl_redraw_all_struct_(iw);
	     swapbuffers();
           } 
	 }
	 /*if t key is pressed translate*/
         else if( getbutton(TKEY) ) {
	   factor = 0.0;
           /*while left button down*/
           while (getbutton(LEFTMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
	     factor = (getvaluator(MOUSEX)-x)*gl.diagonal*0.001;
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     translate( factor,0.0,0.0);
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
  	     gl_redraw_all_struct_(iw);
	     swapbuffers();
           } 
           gl.centre[*iw][0] = gl.centre[*iw][0] - factor;
	 }
	 else { /*rotate*/
           /*while left button down*/
           while (getbutton(LEFTMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     rot(0.3*(getvaluator(MOUSEX)-x), 'x');
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
  	     gl_redraw_all_struct_(iw);
	     swapbuffers();
           } 
	 }
	 getmatrix(gl.modelview[*iw] );
       }
       else if (dev==MIDDLEMOUSE) {
	 x = getvaluator(MOUSEX);
	 /*if t key is pressed translate*/
         if( getbutton(TKEY) ) {
	   factor = 0.0;
           /*while left button down*/
           while (getbutton(MIDDLEMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
	     factor = (getvaluator(MOUSEX)-x)*gl.diagonal*0.001;
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     translate( 0.0,factor,0.0);
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
  	     gl_redraw_all_struct_(iw);
	     swapbuffers();
           } 
           gl.centre[*iw][1] = gl.centre[*iw][1] - factor;
	 }
	 else { /*rotate*/
           /*while left button down*/
	   while (getbutton(MIDDLEMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     rot(0.3*(getvaluator(MOUSEX)-x), 'y');
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
  	     gl_redraw_all_struct_(iw);
	     swapbuffers();
           } 
	 }
	 getmatrix(gl.modelview[*iw] );
       }
       else if (dev==RIGHTMOUSE) {
	 x = getvaluator(MOUSEX);
	 /*if t key is pressed translate*/
         if( getbutton(TKEY) ) {
	   factor = 0.0;
           /*while left button down*/
           while (getbutton(RIGHTMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
	     factor = (getvaluator(MOUSEX)-x)*gl.diagonal*0.001;
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     translate( 0.0,0.0,factor);
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
  	     gl_redraw_all_struct_(iw);
	     swapbuffers();
           } 
           gl.centre[*iw][2] = gl.centre[*iw][2] - factor;
	 }
	 else { /*rotate*/
           /*while left button down*/
	   while (getbutton(RIGHTMOUSE)) {
	     /*poll mouse position and rotate*/
             gl_select_xform_(iw);
             translate( gl.centre[*iw][0], gl.centre[*iw][1],gl.centre[*iw][2]);
	     rot(0.3*(getvaluator(MOUSEX)-x), 'z');
             translate(-gl.centre[*iw][0],-gl.centre[*iw][1],-gl.centre[*iw][2]);
  	     gl_redraw_all_struct_(iw);
	     swapbuffers();
           } 
	 }
	 getmatrix(gl.modelview[*iw] );
       }
       else if (dev==ESCKEY && getbutton(ESCKEY)) cont = 0;
     }   
   } 
   frontbuffer(TRUE);
   depthcue(FALSE);
}


gl_select_xform_(iw)

   /* select input transformation iw */
   int *iw;

{
   loadmatrix(gl.modelview[*iw]);
}


gl_set_animation_flag_(flag)

   int *flag;

{
   animate = *flag;
   if(animate) {
     frame = 0;
     max_segments[frame] = 0;
   }
}


gl_set_animation_frame_(time)

   int *time;

{
   frame = *time;
}


gl_set_colour_rep_(iw,colour_index,red,green,blue)

   /* set colour representation*/
   int *iw;               /* set colour LUT iw */
   int *colour_index;     /* entry of LUT */
   float *red,*green,*blue; /* colour values */

{

   /* changed to RGB mode 9 May 1991 AAY*/
   /* keep copy of colour lut for each ws */
   short r,g,b;
   /*set the appropriate lut number (0..15)*/
   r=(short)(*red * 255.0);
   g=(short)(*green * 255.0);
   b=(short)(*blue * 255.0);
   gl.colourmap[*iw][*colour_index][0] = r;
   gl.colourmap[*iw][*colour_index][1] = g;
   gl.colourmap[*iw][*colour_index][2] = b;
}


gl_set_depth_cue_index_(index)

   /* set index to bundled attributes */
   int *index;

{
   /* attributes are actually set, not bundled, unfortunately */
}


gl_set_depth_cue_rep_(iw,index,type,planes,scales,colour_type,colour)

   /* set depth cue representation */
   int *iw,*index,*type,*colour_type;
   float *planes,*scales,*colour;
{

}


gl_set_edge_index_(index)

   /* set index to bundled attributes */
   int *index;

{
   /* attributes are actually set, not bundled, unfortunately */
}


gl_set_edge_rep_(iw,index,flag,type)

   /* set edge representation for polygons */
   int *iw,*index,*flag,*type;

{
   /* not supported by gl? */
}


gl_set_face_cull_mode_(type)

   /* not supported */
   int *type;

{

}


gl_set_face_dist_mode_(type)

   /* not supported */
   int *type;

{

}


gl_set_fill_area_index_(index)

   /* set index to bundled attributes */
   int *index;

{
   /* attributes are actually set, not bundled, unfortunately */
}


gl_set_fill_rep_(iw,index,pattern,style,colour_index)

   /* set fill representation*/
   int *iw;             /* workstation window */
   int *index;          /* bundled index */
   int *pattern,*style; /* define pattern and style of filled area */
   int *colour_index;   /* colour lut index */

{
   /* a table has to be kept of type constants used in each bundle */
}


gl_set_global_xform_(iw,matrix)

   /* Updates global object transformation */
   int *iw;       /* window or transformation number */
   float *matrix; /* 4x4 matrix */

{

}


gl_set_light_src_state_(n,active,n2,deact)

   /* set light source state */
   int *n,*active,*n2,*deact;

{

}


gl_set_light_src_rep_(iw,index,type,size,data)

   /* set light source representation */
   int *iw,*index,*type,*size,*data;

{

}


gl_set_polyline_index_(iw,index)

   /* set index to bundled attributes */
   int *index,*iw;

{
   int i;
   static short purple[] = {255,0,255};

   /* attributes are actually set, not bundled, unfortunately */
   i = (polyline_rep[*iw-1].colour[*index]);
   if(*index==33) 
     c3s(purple);
   else
     c3s( gl.colourmap[*iw][i] );
   setlinestyle(polyline_rep[*iw-1].type[*index]);
}


gl_set_polyline_rep_(iw,index,type,width,colour_index)

   /* set polyline representation*/
   int *iw;             /* workstation window */
   int *index;          /* bundled index */
   int *type;           /* line type */
   float *width;        /* line width */
   int *colour_index;   /* colour lut index */

{
   /* a table has to be kept of type constants used in each bundle */
   if (*index>MAX_POLYLINES) {
     printf(" too many polyline representations");
     return;
   }
   if (*iw>MAX_WINDOWS) {
     printf(" too many polyline window representations");
     return;
   }
   polyline_rep[*iw-1].type[*index-1]= *type;
   if(*type==3)
     deflinestyle((short)*index,8);
   polyline_rep[*iw-1].width[*index-1]= *width;
   polyline_rep[*iw-1].colour[*index-1]= *colour_index;
}


gl_set_polymarker_index_(iw,index)

   /* set index to bundled attributes */
   int *index,*iw;

{
   int i;
   static short purple[] = {0,255,0};

   /* attributes are actually set, not bundled, unfortunately */
   i = (polymarker_rep[*iw-1].colour[*index]);
   /* c3s( gl.colourmap[*iw][i] ); */
   c3s( purple );
}


gl_set_polymarker_rep_(iw,index,type,size,colour_index)

   /* set polymarker representation*/
   int *iw;             /* workstation window */
   int *index;          /* bundled index */
   int *type;           /* marker type */
   float *size;         /* marker size */
   int *colour_index;   /* colour lut index */

{
   /* a table has to be kept of type constants used in each bundle */
   if (*index>MAX_POLYMARKERS) {
     printf(" too many polymarker representations");
     return;
   }
   if (*iw>MAX_WINDOWS) {
     printf(" too many polymarker window representations");
     return;
   }
   polymarker_rep[*iw-1].type[*index-1]= *type;
   polymarker_rep[*iw-1].size[*index-1]= *size;
   polymarker_rep[*iw-1].colour[*index-1]= *colour_index;
}


gl_set_surface_index_(iw,index)

   /* set index to bundled attributes */
   int *iw,*index;

{
   /* attributes are actually set, not bundled, unfortunately */
   /* bind appropiate lighting source, material and model */
   if (*index == 1 || *index == 2) {
     lmbind(MATERIAL,0);
     lmbind(LIGHT0,0);
     lmbind(LMODEL,0);
   }
   else {
     lmbind(MATERIAL,*index);
     lmbind(LIGHT0,*index);
     lmbind(LIGHT1, 1);
     lmbind(LMODEL,1);
     /*lmbind(MATERIAL, 1);
     lmbind(LIGHT0, 1);
     lmbind(LMODEL, 1);*/
   }
}


gl_set_surface_rep_(iw,index,front_style,back_style,fsi,bsi,
   front_colour_type,front_colour,back_colour_type,back_colour,
   front_shade,back_shade,front_light,back_light,
   f_amb,f_diff,f_spec,fspec_col_type,fspec_col,fspec_exp,f_transparent,
   b_amb,b_diff,b_spec,bspec_col_type,bspec_col,bspec_exp,b_transparent,
   f_approx,f_trim,b_approx,b_trim)

   /* set interior surface representation */
   int *iw,*index,*front_style,*back_style,*fsi,*bsi,*front_colour_type,
       *back_colour_type,*front_shade,*back_shade,*front_light,*back_light,
       *fspec_col_type,*bspec_col_type,*f_approx,*f_trim,*b_approx,*b_trim;
   float *front_colour,*back_colour,
         *f_amb,*f_diff,*f_spec,*fspec_col,*fspec_exp,*f_transparent,
         *b_amb,*b_diff,*b_spec,*bspec_col,*bspec_exp,*b_transparent;

{
  static float lt[] = {
    LCOLOR, 0.6, 0.0, 0.0,
    POSITION, 0.0, 0.1, 0.0, 0.0,
    LMNULL };
  static float mat[] = {
    AMBIENT, 0.7,0.3,0.4,
/*    DIFFUSE, 0.7, 0.5, 0.45, nice flesh brown*/
    DIFFUSE, 0.9, 0.5, 0.6,
    SPECULAR, 0.9,0.5,0.6,
    SHININESS, 30,
    LMNULL };
  static float lm[] = {
    AMBIENT, 0.2,0.2,0.2,
    LOCALVIEWER, 0,
    LMNULL };
  
  /* only use one lighting model at present */
  lmdef(DEFMATERIAL,*index,0,mat);
  lmdef(DEFLIGHT,*index,10,lt);
  lmdef(DEFLMODEL,*index,0,lm); 
/*  lmdef(DEFMATERIAL, 1, 0, NULL);
  lmdef(DEFLIGHT, 1, 0, NULL);*/
  lmdef(DEFLMODEL, 1, 0, NULL);
  lmdef(DEFLIGHT,1,0,NULL);
}


gl_set_text_index_(iw,index)

   /* set index to bundled attributes */
   int *iw;             /* workstation window */
   int *index;

{
   int i;
   /* attributes are actually set, not bundled, unfortunately */
   i = (text_rep[*iw-1].colour[*index]);
   c3s( gl.colourmap[*iw][i] );
}


gl_set_text_geom_(index)

   /* set geometric attributes */
   int *index;

{
   /* attributes are actually set, not bundled, unfortunately */
}


gl_set_text_rep_(iw,index,font,prec,width,spacing,colour_index)

   /* set text representation*/
   int *iw;             /* workstation window */
   int *index;          /* bundled index */
   int *font;           /* text font number */
   float *width;        /* text width */
   float *spacing;      /* text spacing */
   int *colour_index;   /* colour lut index */

{
   /* a table has to be kept of type constants used in each bundle */
   if (*index>MAX_TEXTS) {
     printf(" too many text representations");
     return;
   }
   if (*iw>MAX_WINDOWS) {
     printf(" too many text window representations");
     return;
   }
   text_rep[*iw-1].font[*index-1]= *font;
   text_rep[*iw-1].width[*index-1]= *width;
   text_rep[*iw-1].spacing[*index-1]= *spacing;
   text_rep[*iw-1].colour[*index-1]= *colour_index;
}


gl_set_view_(istatus,iw,mode_proj,npc_clip,angle1,angle2,angle3,
   back_plane_dist,front_plane_dist,fpt,fscale,fshft,proj_ref_point,
   view_plane_dist,view_plane,view_ref_point,view_up,view_port,wind)

   /* update the viewing transform */
   int *istatus;       /*error flag*/
   int *iw;            /*window/transformation number*/
   int *mode_proj;     /*0 for orthographic, 1 for perspective */
   int *npc_clip;      /*not used*/
   float *angle1,*angle2,*angle3; /*view rotation*/
   float *back_plane_dist,*front_plane_dist; /*front and back clipping planes*/
   float *fpt,*fscale,*fshft; /*operates on the projection matrix?*/
   float *proj_ref_point;   /*is the point from which you look*/
   float *view_plane_dist;  /*I forget what this is (usually 0)*/
   float *view_plane;        /*normal to the plane at which you look*/
   float *view_ref_point;    /*the point in the object at which you look*/
   float *view_up;          /*orientation of the viewer ie which way is up*/
   float *view_port;         /*usually not changed from default*/
   float *wind;             /*area of the view plane which gets mapped onto the viewport*/
{
   short a;
   /*combine eval view map, eval view orient and transform view */

   /*first overwrite the projection matrix*/
   if (*mode_proj==0)  /*orthogonal*/
      ortho( wind[0],wind[1],wind[2],wind[3],
	 *front_plane_dist,*back_plane_dist);
   else {
      window( wind[0],wind[1],wind[2],wind[3],
	 *front_plane_dist,*back_plane_dist);
   }

   /*set line of sight*/
   a = atan(view_up[0]/view_up[1]) * 10.0;
   lookat( proj_ref_point[0], proj_ref_point[1], proj_ref_point[2],
      view_ref_point[0], view_ref_point[1], view_ref_point[2], a );
   
   /*add rotate, translate, scale, shift here*/

   getmatrix( gl.modelview[*iw] );
}


gl_set_xform_(iw,fixed_pt,shift,angle1,angle2,angle3,scalefact,istatus)

   /* set window transformation */ 
   int *iw;  /* window and transformation number */
   float *fixed_pt,*shift,*angle1,*angle2,*angle3,*scalefact;
   int *istatus; /* error flag */

{
   /* this must be done after the view is set and not before? */
   rot( *angle1, 'x' );
   rot( *angle2, 'y' );
   rot( *angle3, 'z' );
   translate( shift[0], shift[1], shift[2] );
   scale( scalefact[0], scalefact[1], scalefact[2] );
   getmatrix( gl.modelview[*iw] );
   if(isobj(CURRENT_TRANS)) {
     editobj(CURRENT_TRANS);
     gl_select_xform_(iw);
     closeobj();
   }
}


gl_surface_( iw,shape_flag,facet_data_flag,vertex_data_flag,edge_data_flag,
   egde_vis,colour_type,nvertices,facet_normals,facet_colours,vertices,
   vertex_colours,vertex_normals,nfacets,nvpf,indices)

   /* draw a shaded surface */
   int *iw,*shape_flag,*facet_data_flag,*vertex_data_flag,*edge_data_flag,
       *egde_vis,*colour_type,*nvertices,*facet_colours,*vertex_colours,
       *nfacets,*nvpf,*indices;
   float *facet_normals,*vertices,*vertex_normals;

{
   int nvert,vert,i,j,index;

   /* calculate number of vertices on each side assuming square */
   nvert = sqrt( (double) *nvertices );

   /* if vertex colour data */
   if( *vertex_data_flag == 2 ) {  
     /* do for all strips */
     for( j=1; j<=(nvert-1); j++ ) {
       bgnqstrip();
       /* present vertices in the right order */
       for( i=1; i<=nvert; i++ ) {
         vert = nvert*j + i;
         index = *(vertex_colours+(vert-1));
         c3s( gl.colourmap[*iw][index] );
         v3f( vertices+(vert-1)*3 );
         vert = nvert*(j-1) + i;
         index = *(vertex_colours+(vert-1));
         c3s( gl.colourmap[*iw][index] );
         v3f( vertices+(vert-1)*3 );
       }
       endqstrip();
     }
   }
   else /* use vertex normal data */ {
     /* do for all strips */
     for( j=1; j<=(nvert-1); j++ ) {
       bgnqstrip();
       /* present vertices in the right order */
       for( i=1; i<=nvert; i++ ) {
         vert = nvert*j + i;
         n3f( vertex_normals+(vert-1)*3 );
         v3f( vertices+(vert-1)*3 );
         vert = nvert*(j-1) + i;
         n3f( vertex_normals+(vert-1)*3 );
         v3f( vertices+(vert-1)*3 );
       }
       endqstrip();
     }
   }

}


gl_text_(point,strin,len)

   /* draw text */
   float *point;
   char *strin;
   int *len;

{
   char string[255];
   int i;

   if(*len>99) *len = 99;
   for( i=0;i< *len;i++) 
     string[i] = *(strin+i);
   string[*len]=0;
   cmov(*point,*(point+1),*(point+2));
   charstr(string);
}


gl_visib_(isegnum,ivis,iw)

   /* if ivis=1 the segment is visible, ivis=2 invisible */
   int *isegnum,*ivis,*iw;

{
   if( *ivis==1 )
     iseg[*isegnum]=2;
   else
     iseg[*isegnum]=1;
}


gl_window_(iw,xmin,xmax,ymin,ymax,zmin,zmax)

   /* set world coordinate window */
   int *iw; /* window or transformation number */
   float *xmin,*xmax,*ymin,*ymax,*zmin,*zmax;

{
   gl.wc_window[*iw][0] = *xmin;
   gl.wc_window[*iw][1] = *xmax;
   gl.wc_window[*iw][2] = *ymin;
   gl.wc_window[*iw][3] = *ymax;
   gl.wc_window[*iw][4] = *zmin;
   gl.wc_window[*iw][5] = *zmax;
   gl.wc_window[*iw][6] = sqrt((*xmax-*xmin)*(*xmax-*xmin)+
			   (*ymax-*ymin)*(*ymax-*ymin)+
			   (*zmax-*zmin)*(*zmax-*zmin));
   if(*iw <= 3)
     gl.diagonal = sqrt((*xmax-*xmin)*(*xmax-*xmin)+
      		        (*ymax-*ymin)*(*ymax-*ymin)+
		        (*zmax-*zmin)*(*zmax-*zmin));
   gl.centre[*iw][0] = ( gl.wc_window[*iw][1] + gl.wc_window[*iw][0] )/2.0;
   gl.centre[*iw][1] = ( gl.wc_window[*iw][3] + gl.wc_window[*iw][2] )/2.0;
   gl.centre[*iw][2] = ( gl.wc_window[*iw][5] + gl.wc_window[*iw][4] )/2.0;
}


gl_wkst_viewport_(iw,xmin,xmax,ymin,ymax)

   /*set screen coord viewport of window*/
   int *iw;
   float *xmin,*xmax,*ymin,*ymax;

{
   /* this must be called before the window is opened */
   gl.view_port[*iw][0]=(short)*xmin;
   gl.view_port[*iw][1]=(short)*xmax;
   gl.view_port[*iw][2]=(short)*ymin;
   gl.view_port[*iw][3]=(short)*ymax;
   prefposition( (long)*xmin, (long)*xmax, (long)*ymin, (long)*ymax );

}


gl_wkst_window_(iw,xmin,xmax,ymin,ymax)

   /* set window in eye coords, used for zooming and panning*/
   int *iw; /* window or transformation number */
   float *xmin,*xmax,*ymin,*ymax;

{

}
