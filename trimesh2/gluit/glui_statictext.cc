/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_statictext.cpp - GLUI_StaticText Control


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "glui_stdinc.h"

/****************************** GLUI_StaticText::draw() **********/

void    GLUI_StaticText::draw( int x, int y )
{
  int orig;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();

  draw_text();

  restore_window( orig );
}


/****************************** GLUI_StaticText::set_text() **********/

void    GLUI_StaticText::set_text( const char *text )
{
  int orig = 0;

  if (strcmp(text,name) == 0) return;

  if ( can_draw() )
  {
  /**** Erase old text first *****/
      orig = set_to_glut_window();

  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();
  translate_to_origin();
  erase_text();
  glPopMatrix();
  }

  /* Copy text */
  int old_w = w;
  strncpy((char*)name,text,sizeof(GLUI_String)); 
  update_size();
  if ((w != old_w)&&(glui)) glui->refresh();

  if ( can_draw() )
  {

  /**** Redraw the text in the window ****/
  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();
  translate_to_origin();
  draw_text();
  glPopMatrix();

  restore_window( orig );
}
}


/************************************ GLUI_StaticText::update_size() **********/

void   GLUI_StaticText::update_size( void )
{
  int text_size;

  if ( NOT glui )
    return;

  text_size = string_width( name );

  if ( w < text_size )
    w = text_size;    
}


/****************************** GLUI_StaticText::draw_text() **********/

void    GLUI_StaticText::draw_text( void )
{
  if ( NOT can_draw() )
    return;

  erase_text();
  draw_name( 0, 9 );
}


/****************************** GLUI_StaticText::erase_text() **********/

void    GLUI_StaticText::erase_text( void )
{
  if ( NOT can_draw() )
    return;

  set_to_bkgd_color();
  glDisable( GL_CULL_FACE );
  glBegin( GL_TRIANGLES );
  glVertex2i( 0,0 );   glVertex2i( w, 0 );  glVertex2i( w, h );  
  glVertex2i( 0, 0 );  glVertex2i( w, h );  glVertex2i( 0, h );   
  glEnd();
}



