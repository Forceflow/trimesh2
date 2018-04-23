DESTDIR = ../lib.$(UNAME)
INCLUDES = -I../include -I../include/GL

include ../Makerules

COPTS += -Wno-unused-parameter -Wno-shadow

ifdef DONT_BUILD_FREEGLUT
CFILES =
else
CFILES =	freeglut_callbacks.c \
		freeglut_cursor.c \
		freeglut_display.c \
		freeglut_ext.c \
		freeglut_font.c \
		freeglut_font_data.c \
		freeglut_gamemode.c \
		freeglut_geometry.c \
		freeglut_glutfont_definitions.c \
		freeglut_init.c \
		freeglut_input_devices.c \
		freeglut_joystick.c \
		freeglut_main.c \
		freeglut_menu.c \
		freeglut_misc.c \
		freeglut_overlay.c \
		freeglut_spaceball.c \
		freeglut_state.c \
		freeglut_stroke_mono_roman.c \
		freeglut_stroke_roman.c \
		freeglut_structure.c \
		freeglut_teapot.c \
		freeglut_videoresize.c \
		freeglut_window.c \
		freeglut_xinput.c
endif

CCFILES =	glui.cc \
		glui_add_controls.cc \
		glui_algebra3.cc \
		glui_arcball.cc \
		glui_bitmap_img_data.cc \
		glui_bitmaps.cc \
		glui_button.cc \
		glui_checkbox.cc \
		glui_column.cc \
		glui_control.cc \
		glui_edittext.cc \
		glui_listbox.cc \
		glui_mouse_iaction.cc \
		glui_node.cc \
		glui_panel.cc \
		glui_quaternion.cc \
		glui_radio.cc \
		glui_rollout.cc \
		glui_rotation.cc \
		glui_separator.cc \
		glui_slider.cc \
		glui_spinner.cc \
		glui_statictext.cc \
		glui_translation.cc

OFILES = $(addprefix $(OBJDIR)/,$(CFILES:.c=.o) $(CCFILES:.cc=.o))
LIB = $(DESTDIR)/libgluit.a
default: $(LIB)

$(LIB): $(OFILES)
	$(STATICLIB)

clean:
	-rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d
	-rmdir $(OBJDIR)

spotless: clean
	-rm -f $(LIB)
	-rmdir $(DESTDIR)

