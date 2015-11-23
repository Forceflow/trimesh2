/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_slider - GLUI_Slider control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "glui_stdinc.h"

#define GLUI_SLIDER_FONT_HEIGHT					9
#define GLUI_SLIDER_FONT_DROP						3
#define GLUI_SLIDER_FONT_FULL_HEIGHT			(GLUI_SLIDER_FONT_HEIGHT + GLUI_SLIDER_FONT_DROP)
#define GLUI_SLIDER_FONT_MID_HEIGHT				4

#define GLUI_SLIDER_NAME_INDENT					6
#define GLUI_SLIDER_NAME_SIDE_BORDER			2
#define GLUI_SLIDER_NAME_TOP_BORDER				2
#define GLUI_SLIDER_NAME_BOTTOM_BORDER			0

#define GLUI_SLIDER_VAL_TOP_BORDER				3
#define GLUI_SLIDER_VAL_SIDE_BORDER				4

#define GLUI_SLIDER_KNOB_HALF_WIDTH				3
#define GLUI_SLIDER_KNOB_HALF_HEIGHT			8
#define GLUI_SLIDER_KNOB_SIDE_BORDER			4
#define GLUI_SLIDER_KNOB_TOP_BORDER				2
#define GLUI_SLIDER_KNOB_BOTTOM_BORDER			2
#define GLUI_SLIDER_KNOB_WHITE					254
#define GLUI_SLIDER_KNOB_L_GREY					230
#define GLUI_SLIDER_KNOB_GREY						218
#define GLUI_SLIDER_KNOB_D_GREY					115
#define GLUI_SLIDER_KNOB_BLACK					0


#define GLUI_SLIDER_MAX_VAL_STRING_SIZE      200
#define GLUI_SLIDER_TRUNCATE_STRING				"~"


/********************** GLUI_Slider::mouse_down_handler() ******/

int    GLUI_Slider::mouse_down_handler( int local_x, int local_y )
{
	pressed = true;

	//printf("%d %d\n",local_x,local_y);

	update_val(local_x - x_abs, local_y - y_abs);
	draw_translated_active_area();
	if( glui ) glui->post_update_main_gfx();
	execute_callback();

	return false;
}


/**************************** GLUI_Slider::mouse_up_handler() */

int    GLUI_Slider::mouse_up_handler( int local_x, int local_y, int inside )
{
	pressed = false;

	//printf("%d %d %d\n",local_x,local_y,inside);

	update_val(local_x - x_abs, local_y - y_abs);
	draw_translated_active_area();
	if( glui ) glui->post_update_main_gfx();
	execute_callback();

	return false;
}


/****************************** GLUI_Slider::mouse_held_down_handler() ******/

int    GLUI_Slider::mouse_held_down_handler( int local_x, int local_y, int inside)
{  

	//printf("%d %d %d\n",local_x,local_y,inside);

	update_val(local_x - x_abs, local_y - y_abs);
	draw_translated_active_area();
	if( glui ) glui->post_update_main_gfx();
	execute_callback();

	return false;
}

/****************************** GLUI_Slider::special_handler() **********/

int    GLUI_Slider::special_handler( int key,int mods )
{
	int grads;

	int min_x, max_x, wid_x;

	int g;
	float g_f;

	int new_int;

	min_x = 2 + GLUI_SLIDER_KNOB_SIDE_BORDER + GLUI_SLIDER_KNOB_HALF_WIDTH;
	max_x = w - 2 - GLUI_SLIDER_KNOB_HALF_WIDTH - GLUI_SLIDER_KNOB_SIDE_BORDER;
	wid_x = max_x - min_x + 1;

	//The arrows adjust the slider's value.
	//Without mods (shift, ctrl, alt), the current
	//value is rounded to the nearest grad and then
	//one grad is added/subtracted.  However, shift,
	//ctrl, alt modify the grad size by 1/2,1/10, 
	//1/100. If this is an int slider, the min
	//change is 1.
	//Note, by default, grads=0 which takes the
	//num of grads to be equal to the pixel
	//width of the slider...

	if (graduations == 0)
		grads = wid_x;
	else
		grads = graduations;

	if (mods == 0)
	{
		//Use unmodified grads
	}
	else if (mods & GLUT_ACTIVE_SHIFT)
	{
		grads = 2*(grads-1) + 1;
	}
	else if (mods & GLUT_ACTIVE_CTRL)
	{
		grads = 10*(grads-1) + 1;
	}
	else if (mods & GLUT_ACTIVE_ALT)
	{
		grads = 100*(grads-1) + 1;
	}
	else
	{
		return false;
	}

	switch (data_type)
	{
		case GLUI_SLIDER_INT:
			if (int_low != int_high)
			{
				if (int_val == int_high)
					g = grads-1;
				else
					g = ((int)(((double)grads)*((double)(int_val-int_low))/((double)(int_high - int_low))));
			}
			else
				g = 0;
		break;
		case GLUI_SLIDER_FLOAT:
			if (float_low != float_high)
			{
				if (float_val == float_high)
					g = grads-1;
				else
					g = ((int)(((double)grads)*((double)(float_val-float_low))/((double)(float_high - float_low))));
			}
			else
				g = 0;
		break;
		default:
			fprintf(stderr,"GLUI_Slider::upate_knob - Impossible data type!\n");
			abort();
	}

	switch (key)
	{
		case GLUT_KEY_RIGHT:
			g += 1;
			if (g > (grads-1))
				g = grads-1;
		break;
		case GLUT_KEY_LEFT:
			g -= 1;
			if (g < 0)
				g = 0;
		break;
		default:
			return false;
		break;
	}

	g_f = (float) g / (grads - 1);
	
	switch (data_type)
	{
		case GLUI_SLIDER_INT:
			new_int = int(int_low + g_f * (int_high - int_low));

			if (new_int == int_val)
			{
				switch (key)
				{
					case GLUT_KEY_RIGHT:
						new_int += 1;
						if (new_int > int_high)
							new_int = int_high;
					break;
					case GLUT_KEY_LEFT:
						new_int -= 1;
						if (new_int < int_low)
							new_int = int_low;
					break;
				}
			}
			set_int_val(new_int);     
			set_float_val((float)(int_val));
		break;
		case GLUI_SLIDER_FLOAT:
			set_float_val(float_low + g_f*(float_high-float_low));
			set_int_val(int(float_val));
		break;
		default:
			fprintf(stderr,"GLUI_Slider::upate_knob - Impossible data type!\n");
			abort();
	}

	draw_translated_active_area();
	if( glui ) glui->post_update_main_gfx();
	execute_callback();

	return false;
}

/****************************** GLUI_Slider::update_knob() **********/

void    GLUI_Slider::update_val( int x, int y )
{

	int grads;

	int min_x, max_x, wid_x;

	int g;
	double g_f;

	min_x = 2 + GLUI_SLIDER_KNOB_SIDE_BORDER + GLUI_SLIDER_KNOB_HALF_WIDTH;
	max_x = w - 2 - GLUI_SLIDER_KNOB_HALF_WIDTH - GLUI_SLIDER_KNOB_SIDE_BORDER;
	wid_x = max_x - min_x + 1;

	if (graduations == 0)
		grads = wid_x;
	else
		grads = graduations;

	if (x < min_x)
		x = min_x;
	else if (x > max_x)
		x = max_x;

	if (x == max_x)
		g = grads-1;
	else
		g = (int) ((double)(((double)grads)*((double)(x-min_x))/((double)(max_x - min_x))));

	g_f = ((double)g)/((double)(grads-1));
	
	switch (data_type)
	{
		case GLUI_SLIDER_INT:
			set_int_val((int)(((double)int_low)+0.5+((double)(g_f*((double)(int_high-int_low))))));     
			set_float_val((float)(int_val));
		break;
		case GLUI_SLIDER_FLOAT:
			set_float_val((float)(((double)float_low) + g_f*((double)(float_high-float_low))));
			set_int_val((int)(float_val));
		break;
		default:
			fprintf(stderr,"GLUI_Slider::upate_knob - Impossible data type!\n");
			abort();
	}

}
 
/****************************** GLUI_Slider::draw() **********/

void    GLUI_Slider::draw( int x, int y )
{
  int orig;

  if ( NOT glui )
    return;

	orig = set_to_glut_window();


	draw_emboss_box(	0, 
							w,
							GLUI_SLIDER_NAME_TOP_BORDER +
							GLUI_SLIDER_FONT_HEIGHT - 1 - 
							GLUI_SLIDER_FONT_MID_HEIGHT, 
							h );

	draw_bkgd_box( GLUI_SLIDER_NAME_INDENT-1, 
						GLUI_SLIDER_NAME_INDENT + 
						string_width(name.string) + 
						2*GLUI_SLIDER_NAME_SIDE_BORDER - 1, 
						0, 
						0 + 
						GLUI_SLIDER_FONT_FULL_HEIGHT + 
						GLUI_SLIDER_NAME_TOP_BORDER + 
						GLUI_SLIDER_NAME_BOTTOM_BORDER);

	draw_name(	GLUI_SLIDER_NAME_INDENT + 
					GLUI_SLIDER_NAME_SIDE_BORDER, 
					GLUI_SLIDER_FONT_HEIGHT-1 + 
					GLUI_SLIDER_NAME_TOP_BORDER);

	draw_active_area();

	restore_window(orig);

	draw_active_box(	GLUI_SLIDER_NAME_INDENT, 
							GLUI_SLIDER_NAME_INDENT +
							string_width(name.string) + 
							2*GLUI_SLIDER_NAME_SIDE_BORDER - 1,
							0, 
							GLUI_SLIDER_FONT_FULL_HEIGHT + 
							GLUI_SLIDER_NAME_TOP_BORDER + 
							GLUI_SLIDER_NAME_BOTTOM_BORDER - 1);

}


/************************************ GLUI_Slider::update_size() **********/

void   GLUI_Slider::update_size( void )
{
	int min_w;
	
	if ( NOT glui )
		return;

	h =	GLUI_SLIDER_NAME_TOP_BORDER +
			GLUI_SLIDER_FONT_HEIGHT - 
			GLUI_SLIDER_FONT_MID_HEIGHT +
			GLUI_SLIDER_VAL_TOP_BORDER + 
			GLUI_SLIDER_FONT_HEIGHT  +
			GLUI_SLIDER_KNOB_TOP_BORDER +
			GLUI_SLIDER_KNOB_HALF_HEIGHT +
			1 +
			GLUI_SLIDER_KNOB_HALF_HEIGHT +
			GLUI_SLIDER_KNOB_BOTTOM_BORDER +
			2;

	min_w =	GLUI_SLIDER_NAME_INDENT +
				string_width(name.string) + 
				2*GLUI_SLIDER_NAME_SIDE_BORDER +
				string_width(GLUI_SLIDER_TRUNCATE_STRING) + 
				string_width("0000000") +
				2*GLUI_SLIDER_VAL_SIDE_BORDER + 
				2 - 1; //+2 for border width at right
				

	if (w<min_w)
		w = min_w;

}


/****************************** GLUI_Mouse_Interaction::draw_active_area_translated() **********/

void    GLUI_Slider::draw_translated_active_area( void )
{
	int orig;
	int win_h , win_w;

	if ( NOT glui )
		return;

	orig = set_to_glut_window();
	win_h = glutGet( GLUT_WINDOW_HEIGHT );
	win_w = glutGet(GLUT_WINDOW_WIDTH);
	
	glMatrixMode( GL_MODELVIEW );
	glPushMatrix();
	glLoadIdentity();
	glTranslatef( (float) win_w/2.0f, (float) win_h/2.0f, 0.0f );
	glRotatef( 180.0f, 0.0f, 1.0f, 0.0f );
	glRotatef( 180.0f, 0.0f, 0.0f, 1.0f );
	glTranslatef( (float) -win_w/2.0f, (float) -win_h/2.0f, 0.0f );
	glTranslatef( (float) this->x_abs + .5f, (float) this->y_abs + .5f, 0.0f );

	draw_active_area();
	
	glMatrixMode( GL_MODELVIEW );
	glPopMatrix();

	restore_window(orig);
}

/****************************** GLUI_Mouse_Interaction::draw_active_area() **********/

void    GLUI_Slider::draw_active_area( void )
{
	draw_val();
	draw_slider();
}


/************************************ GLUI_Slider::draw_val() **********/

void   GLUI_Slider::draw_val( void )
{
	int max_w, i;

	char buf[GLUI_SLIDER_MAX_VAL_STRING_SIZE];

	switch (data_type)
	{
		case GLUI_SLIDER_INT:
			sprintf(buf,"%d", int_val);
		break;
		case GLUI_SLIDER_FLOAT:
			sprintf(buf,"%.3g", float_val);
		break;
		default:
			fprintf(stderr,"GLUI_Slider::draw_val - Impossible data type!\n");
			abort();
	}

	max_w =	w +1 -
				(GLUI_SLIDER_NAME_INDENT +
				string_width(name.string) + 
				2*GLUI_SLIDER_NAME_SIDE_BORDER +
				2*GLUI_SLIDER_VAL_SIDE_BORDER + 2);

	if (max_w < string_width(GLUI_SLIDER_TRUNCATE_STRING))
	{
		fprintf(stderr,"GLUI_Slider::draw_val - Impossible max_w!!!!\n");
		abort();	
	}

	if (string_width(buf) > max_w)
	{

		max_w -= string_width(GLUI_SLIDER_TRUNCATE_STRING);

		i = strlen(buf)-1;
		while ((string_width(buf) > max_w) && (i>=0))
		{
			buf[i--] = 0;
		}

		sprintf(&buf[i+1],"%s",GLUI_SLIDER_TRUNCATE_STRING);
	}

   draw_bkgd_box(	
						GLUI_SLIDER_NAME_INDENT +
						string_width(name.string) + 
						2*GLUI_SLIDER_NAME_SIDE_BORDER +
						GLUI_SLIDER_VAL_SIDE_BORDER - 1,

						w-1 - GLUI_SLIDER_VAL_SIDE_BORDER,

						GLUI_SLIDER_NAME_TOP_BORDER +
						GLUI_SLIDER_FONT_HEIGHT - 
						GLUI_SLIDER_FONT_MID_HEIGHT +
						GLUI_SLIDER_VAL_TOP_BORDER + 
						GLUI_SLIDER_FONT_HEIGHT + 1, 

						GLUI_SLIDER_NAME_TOP_BORDER +
						GLUI_SLIDER_FONT_HEIGHT - 
						GLUI_SLIDER_FONT_MID_HEIGHT +
						GLUI_SLIDER_VAL_TOP_BORDER + 1 );

	glColor3b( 0, 0, 0 );
   glRasterPos2i(	w - 
						string_width(buf) - 
						GLUI_SLIDER_VAL_SIDE_BORDER,
						GLUI_SLIDER_NAME_TOP_BORDER +
						GLUI_SLIDER_FONT_HEIGHT - 
						GLUI_SLIDER_FONT_MID_HEIGHT +
						GLUI_SLIDER_VAL_TOP_BORDER + 
						GLUI_SLIDER_FONT_HEIGHT); 
   draw_string(buf);

}

/************************************ GLUI_Slider::draw_slider() **********/

void   GLUI_Slider::draw_slider( void )
{

	int min_x;
	int max_x;
	int wid_x;

	float min_val;
	float max_val;
	float val;

	int x_pos;

   draw_bkgd_box(	
						2 + GLUI_SLIDER_KNOB_SIDE_BORDER - 1,
						w-2 - GLUI_SLIDER_KNOB_SIDE_BORDER,
						GLUI_SLIDER_NAME_TOP_BORDER +
						GLUI_SLIDER_FONT_HEIGHT - 
						GLUI_SLIDER_FONT_MID_HEIGHT +
						GLUI_SLIDER_VAL_TOP_BORDER + 
						GLUI_SLIDER_FONT_HEIGHT + 
						GLUI_SLIDER_KNOB_TOP_BORDER+1,
						GLUI_SLIDER_NAME_TOP_BORDER +
						GLUI_SLIDER_FONT_HEIGHT - 
						GLUI_SLIDER_FONT_MID_HEIGHT +
						GLUI_SLIDER_VAL_TOP_BORDER + 
						GLUI_SLIDER_FONT_HEIGHT + 
						GLUI_SLIDER_KNOB_TOP_BORDER +
						1+ 2*GLUI_SLIDER_KNOB_HALF_HEIGHT + 1);


	draw_emboss_box(				2 +
										GLUI_SLIDER_KNOB_SIDE_BORDER +
										GLUI_SLIDER_KNOB_HALF_WIDTH,
										w - 2 -
										GLUI_SLIDER_KNOB_HALF_WIDTH -
										GLUI_SLIDER_KNOB_SIDE_BORDER,
										GLUI_SLIDER_NAME_TOP_BORDER +
										GLUI_SLIDER_FONT_HEIGHT - 
										GLUI_SLIDER_FONT_MID_HEIGHT +
										GLUI_SLIDER_VAL_TOP_BORDER + 
										GLUI_SLIDER_FONT_HEIGHT  +
										GLUI_SLIDER_KNOB_TOP_BORDER +
										GLUI_SLIDER_KNOB_HALF_HEIGHT,
										GLUI_SLIDER_NAME_TOP_BORDER +
										GLUI_SLIDER_FONT_HEIGHT - 
										GLUI_SLIDER_FONT_MID_HEIGHT +
										GLUI_SLIDER_VAL_TOP_BORDER + 
										GLUI_SLIDER_FONT_HEIGHT + 
										GLUI_SLIDER_KNOB_TOP_BORDER +
										GLUI_SLIDER_KNOB_HALF_HEIGHT + 3);


	min_x = 2 + GLUI_SLIDER_KNOB_SIDE_BORDER + GLUI_SLIDER_KNOB_HALF_WIDTH;
	max_x = w - 2 - GLUI_SLIDER_KNOB_HALF_WIDTH - GLUI_SLIDER_KNOB_SIDE_BORDER;
	wid_x = max_x - min_x + 1;

	switch (data_type)
	{
		case GLUI_SLIDER_INT:
			min_val = (float)int_low;
			max_val = (float)int_high;
			val = (float)int_val;
		break;
		case GLUI_SLIDER_FLOAT:
			min_val = float_low;
			max_val = float_high;
			val = float_val;
		break;
		default:
			fprintf(stderr,"GLUI_Slider::draw_val - Impossible data type!\n");
			abort();
		break;
	}

	if (max_val == min_val)
		x_pos = min_x;
	else
	{
		if (val == max_val)
			x_pos = max_x;
		else
			x_pos = (int) (min_x + wid_x*((val-min_val)/(max_val-min_val)));
	}

	draw_knob(x_pos,0,0,0,0,glui->bkgd_color.r, glui->bkgd_color.g, glui->bkgd_color.b);
  	draw_knob(x_pos,1,0,1,0, GLUI_SLIDER_KNOB_BLACK, GLUI_SLIDER_KNOB_BLACK, GLUI_SLIDER_KNOB_BLACK);
  	draw_knob(x_pos,0,1,0,1, GLUI_SLIDER_KNOB_WHITE, GLUI_SLIDER_KNOB_WHITE, GLUI_SLIDER_KNOB_WHITE);
	if (pressed)
		draw_knob(x_pos,1,1,1,1, GLUI_SLIDER_KNOB_GREY, GLUI_SLIDER_KNOB_GREY, GLUI_SLIDER_KNOB_GREY);
	else
		draw_knob(x_pos,1,1,1,1,glui->bkgd_color.r, glui->bkgd_color.g, glui->bkgd_color.b);
  	draw_knob(x_pos,2,1,2,1, GLUI_SLIDER_KNOB_D_GREY, GLUI_SLIDER_KNOB_D_GREY, GLUI_SLIDER_KNOB_D_GREY);
  	if (pressed)
		draw_knob(x_pos,2,2,2,2, GLUI_SLIDER_KNOB_L_GREY, GLUI_SLIDER_KNOB_L_GREY, GLUI_SLIDER_KNOB_L_GREY);	
	else
		draw_knob(x_pos,2,2,2,2, GLUI_SLIDER_KNOB_GREY, GLUI_SLIDER_KNOB_GREY, GLUI_SLIDER_KNOB_GREY);
  	

}
/******************************** GLUI_Slider::draw_knob() **********/

void		GLUI_Slider::draw_knob( 
				int x, 
				int off_l, 
				int off_r, 
				int off_t, 
				int off_b, 
				unsigned char r,
				unsigned char g,
				unsigned char b)
{

	float col_r = ((float)r)/255.0f;
	float col_g = ((float)g)/255.0f;
	float col_b = ((float)b)/255.0f;

	draw_box(													
		x - GLUI_SLIDER_KNOB_HALF_WIDTH - 1 +		
		off_l,													
		x + GLUI_SLIDER_KNOB_HALF_WIDTH -			
		off_r,													
		GLUI_SLIDER_NAME_TOP_BORDER +						
		GLUI_SLIDER_FONT_HEIGHT -							
		GLUI_SLIDER_FONT_MID_HEIGHT +						
		GLUI_SLIDER_VAL_TOP_BORDER +						
		GLUI_SLIDER_FONT_HEIGHT  +							
		GLUI_SLIDER_KNOB_TOP_BORDER + 1 +				
		off_t,													
		GLUI_SLIDER_NAME_TOP_BORDER +						
		GLUI_SLIDER_FONT_HEIGHT -							
		GLUI_SLIDER_FONT_MID_HEIGHT +						
		GLUI_SLIDER_VAL_TOP_BORDER +						
		GLUI_SLIDER_FONT_HEIGHT  +							
		GLUI_SLIDER_KNOB_TOP_BORDER +						
		GLUI_SLIDER_KNOB_HALF_HEIGHT + 1 +				
		GLUI_SLIDER_KNOB_HALF_HEIGHT + 1 -				
		off_b,													
		col_r,
		col_g,
		col_b);


}

/********************************* GLUI_Slider::set_float_limits() *********/

void GLUI_Slider::set_float_limits( float low, float high)
{

	if (low > high)
	{
		fprintf(stderr,"GLUI_Slider::set_float_limits - Ignoring: low > hi\n");
		return;
	}

	float_low   = low;
	float_high  = high;
  
	if ( NOT IN_BOUNDS( float_val, float_low, float_high ))
	{
		set_float_val( float_low );
		set_int_val( (int)float_val );
	}

	int_low     = (int) float_low;
	int_high    = (int) float_high;
}


/*********************************** GLUI_Slider::set_int_limits() *********/

void   GLUI_Slider::set_int_limits( int low, int high)
{
	if (low > high)
	{
		fprintf(stderr,"GLUI_Slider::set_int_limits - Ignoring: low > hi\n");
		return;
	}
	
	int_low     = low;
	int_high    = high;

	if ( NOT IN_BOUNDS( int_val, int_low, int_high ))
	{
		set_int_val( int_low );
		set_float_val( (float)int_val );
	}

	float_low   = (float) int_low;
	float_high  = (float) int_high;
}


/*********************************** GLUI_Slider::set_num_graduations() *********/

void   GLUI_Slider::set_num_graduations( int g )
{

	if (g >= 0)
		graduations = g;

}

/************** GLUI_Slider::GLUI_Slider() ********************/

GLUI_Slider::GLUI_Slider( void ) 
{
	sprintf( name, "Slider: %p", this );

	type							= GLUI_CONTROL_SLIDER;
	w								= GLUI_SLIDER_WIDTH;
	h								= GLUI_SLIDER_HEIGHT;
	can_activate				= true;
	live_type					= GLUI_LIVE_NONE;
	alignment					= GLUI_ALIGN_CENTER;


	int_low						= GLUI_SLIDER_LO;
	int_high						= GLUI_SLIDER_HI;
	float_low					= GLUI_SLIDER_LO;
	float_high					= GLUI_SLIDER_HI;

	pressed						= false;

	graduations					= 0;
}
