/*
Szymon Rusinkiewicz
Princeton University

GLManager.cc
Manager for modern OpenGL buffers and shaders.
*/

#include "GLManager.h"

#include <cstddef>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <cctype>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>

#if defined(_WIN32)
# define WIN32_LEAN_AND_MEAN
# include <windows.h>
# include <GL/gl.h>
#elif defined(__APPLE__)
# include <OpenGL/gl.h>
# include <dlfcn.h>
#else
# include <GL/gl.h>
# include <GL/glx.h>
#endif

using namespace std;

namespace trimesh {

// General #defines
#ifndef APIENTRY
# if defined(GLAPIENTRY)
#  define APIENTRY GLAPIENTRY
# elif defined(GL_APIENTRY)
#  define APIENTRY GL_APIENTRY
# elif defined(_WIN32)
#  define APIENTRY __stdcall
# else
#  define APIENTRY
# endif
#endif
#ifndef APIENTRYP
# define APIENTRYP APIENTRY *
#endif
#define GL_VBO_OFFSET(i) ((GLvoid *)(size_t)(i))


// OpenGL defines beyond GL_VERSION_1_1

// GL_VERSION_1_2
#define GL_UNSIGNED_BYTE_3_3_2            0x8032
#define GL_UNSIGNED_SHORT_4_4_4_4         0x8033
#define GL_UNSIGNED_SHORT_5_5_5_1         0x8034
#define GL_UNSIGNED_INT_8_8_8_8           0x8035
#define GL_UNSIGNED_INT_10_10_10_2        0x8036
#define GL_BGR                            0x80E0
#define GL_BGRA                           0x80E1
#define GL_UNSIGNED_BYTE_2_3_3_REV        0x8362
#define GL_UNSIGNED_SHORT_5_6_5           0x8363
#define GL_UNSIGNED_SHORT_5_6_5_REV       0x8364
#define GL_UNSIGNED_SHORT_4_4_4_4_REV     0x8365
#define GL_UNSIGNED_SHORT_1_5_5_5_REV     0x8366
#define GL_UNSIGNED_INT_8_8_8_8_REV       0x8367
#define GL_UNSIGNED_INT_2_10_10_10_REV    0x8368
// GL_VERSION_1_3
#define GL_TEXTURE0                       0x84C0
#define GL_ACTIVE_TEXTURE                 0x84E0
// GL_VERSION_1_5
#define GL_ARRAY_BUFFER                   0x8892
#define GL_ELEMENT_ARRAY_BUFFER           0x8893
#define GL_ELEMENT_ARRAY_BUFFER_BINDING   0x8895
#define GL_STATIC_DRAW                    0x88E4
// GL_VERSION_2_0
#define GL_FRAGMENT_SHADER                0x8B30
#define GL_VERTEX_SHADER                  0x8B31
#define GL_FLOAT_VEC2                     0x8B50
#define GL_FLOAT_VEC3                     0x8B51
#define GL_FLOAT_VEC4                     0x8B52
#define GL_INT_VEC2                       0x8B53
#define GL_INT_VEC3                       0x8B54
#define GL_INT_VEC4                       0x8B55
#define GL_FLOAT_MAT2                     0x8B5A
#define GL_FLOAT_MAT3                     0x8B5B
#define GL_FLOAT_MAT4                     0x8B5C
#define GL_SAMPLER_2D                     0x8B5E
#define GL_COMPILE_STATUS                 0x8B81
#define GL_LINK_STATUS                    0x8B82
#define GL_VALIDATE_STATUS                0x8B83
#define GL_INFO_LOG_LENGTH                0x8B84
#define GL_ACTIVE_UNIFORMS                0x8B86
#define GL_ACTIVE_UNIFORM_MAX_LENGTH      0x8B87
#define GL_ACTIVE_ATTRIBUTES              0x8B89
#define GL_ACTIVE_ATTRIBUTE_MAX_LENGTH    0x8B8A
// GL_VERSION_3_0
#define GL_HALF_FLOAT                     0x140B
#define GL_RG                             0x8227
#define GL_R16F                           0x822D
#define GL_R32F                           0x822E
#define GL_RG16F                          0x822F
#define GL_RG32F                          0x8230
#define GL_DEPTH_STENCIL                  0x84F9
#define GL_RGBA32F                        0x8814
#define GL_RGB32F                         0x8815
#define GL_RGBA16F                        0x881A
#define GL_RGB16F                         0x881B
#define GL_CLAMP_VERTEX_COLOR             0x891A
#define GL_CLAMP_FRAGMENT_COLOR           0x891B
#define GL_CLAMP_READ_COLOR               0x891C
#define GL_FRAMEBUFFER_BINDING            0x8CA6
#define GL_COLOR_ATTACHMENT0              0x8CE0
#define GL_DEPTH_ATTACHMENT               0x8D00
#define GL_FRAMEBUFFER                    0x8D40
// GL_VERSION_3_1
#define GL_PRIMITIVE_RESTART              0x8F9D
// GL_VERSION_4_3
#define GL_PRIMITIVE_RESTART_FIXED_INDEX  0x8D69
// GL_NV_primitive_restart
#define GL_PRIMITIVE_RESTART_NV           0x8558


// OpenGL typedefs beyond GL_VERSION_1_1

// GL_VERSION_1_5
typedef ptrdiff_t GLsizeiptr;
typedef ptrdiff_t GLintptr;
// GL_VERSION_2_0
typedef char GLchar;


// OpenGL functions beyond GL_VERSION_1_1 on Win32, or GL_VERSION_1_3 on others

// GL_VERSION_1_3
#ifdef _WIN32
typedef void (APIENTRYP glActiveTexture_t) (GLenum texture);
static glActiveTexture_t glActiveTexture;
#endif
// GL_VERSION_1_5
typedef void (APIENTRYP glGenBuffers_t) (GLsizei n, GLuint *buffers);
static glGenBuffers_t glGenBuffers;
typedef void (APIENTRYP glBufferData_t) (GLenum target, GLsizeiptr size, const void *data, GLenum usage);
static glBufferData_t glBufferData;
typedef void (APIENTRYP glBindBuffer_t) (GLenum target, GLuint buffer);
static glBindBuffer_t glBindBuffer;
typedef void (APIENTRYP glDeleteBuffers_t) (GLsizei n, const GLuint *buffers);
static glDeleteBuffers_t glDeleteBuffers;
// GL_VERSION_2_0
typedef GLuint (APIENTRYP glCreateShader_t) (GLenum type);
static glCreateShader_t glCreateShader;
typedef void (APIENTRYP glShaderSource_t) (GLuint shader, GLsizei count, const GLchar *const*string, const GLint *length);
static glShaderSource_t glShaderSource;
typedef void (APIENTRYP glCompileShader_t) (GLuint shader);
static glCompileShader_t glCompileShader;
typedef void (APIENTRYP glDeleteShader_t) (GLuint shader);
static glDeleteShader_t glDeleteShader;
typedef void (APIENTRYP glGetShaderiv_t) (GLuint shader, GLenum pname, GLint *params);
static glGetShaderiv_t glGetShaderiv;
typedef void (APIENTRYP glGetShaderInfoLog_t) (GLuint shader, GLsizei bufSize, GLsizei *length, GLchar *infoLog);
static glGetShaderInfoLog_t glGetShaderInfoLog;
typedef GLuint (APIENTRYP glCreateProgram_t) (void);
static glCreateProgram_t glCreateProgram;
typedef void (APIENTRYP glAttachShader_t) (GLuint program, GLuint shader);
static glAttachShader_t glAttachShader;
typedef void (APIENTRYP glLinkProgram_t) (GLuint program);
static glLinkProgram_t glLinkProgram;
typedef void (APIENTRYP glValidateProgram_t) (GLuint program);
static glValidateProgram_t glValidateProgram;
typedef void (APIENTRYP glUseProgram_t) (GLuint program);
static glUseProgram_t glUseProgram;
typedef void (APIENTRYP glDeleteProgram_t) (GLuint program);
static glDeleteProgram_t glDeleteProgram;
typedef void (APIENTRYP glGetProgramiv_t) (GLuint program, GLenum pname, GLint *params);
static glGetProgramiv_t glGetProgramiv;
typedef void (APIENTRYP glGetProgramInfoLog_t) (GLuint program, GLsizei bufSize, GLsizei *length, GLchar *infoLog);
static glGetProgramInfoLog_t glGetProgramInfoLog;
typedef void (APIENTRYP glUniform1f_t) (GLint location, GLfloat v0);
static glUniform1f_t glUniform1f;
typedef void (APIENTRYP glUniform2f_t) (GLint location, GLfloat v0, GLfloat v1);
static glUniform2f_t glUniform2f;
typedef void (APIENTRYP glUniform3f_t) (GLint location, GLfloat v0, GLfloat v1, GLfloat v2);
static glUniform3f_t glUniform3f;
typedef void (APIENTRYP glUniform4f_t) (GLint location, GLfloat v0, GLfloat v1, GLfloat v2, GLfloat v3);
static glUniform4f_t glUniform4f;
typedef void (APIENTRYP glUniform1i_t) (GLint location, GLint v0);
static glUniform1i_t glUniform1i;
typedef void (APIENTRYP glUniform2i_t) (GLint location, GLint v0, GLint v1);
static glUniform2i_t glUniform2i;
typedef void (APIENTRYP glUniform3i_t) (GLint location, GLint v0, GLint v1, GLint v2);
static glUniform3i_t glUniform3i;
typedef void (APIENTRYP glUniform4i_t) (GLint location, GLint v0, GLint v1, GLint v2, GLint v3);
static glUniform4i_t glUniform4i;
typedef void (APIENTRYP glUniformMatrix2fv_t) (GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
static glUniformMatrix2fv_t glUniformMatrix2fv;
typedef void (APIENTRYP glUniformMatrix3fv_t) (GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
static glUniformMatrix3fv_t glUniformMatrix3fv;
typedef void (APIENTRYP glUniformMatrix4fv_t) (GLint location, GLsizei count, GLboolean transpose, const GLfloat *value);
static glUniformMatrix4fv_t glUniformMatrix4fv;
typedef void (APIENTRYP glGetActiveUniform_t) (GLuint program, GLuint index, GLsizei bufSize, GLsizei *length, GLint *size, GLenum *type, GLchar *name);
static glGetActiveUniform_t glGetActiveUniform;
typedef GLint (APIENTRYP glGetUniformLocation_t) (GLuint program, const GLchar *name);
static glGetUniformLocation_t glGetUniformLocation;
typedef void (APIENTRYP glVertexAttrib1f_t) (GLuint index, GLfloat x);
static glVertexAttrib1f_t glVertexAttrib1f;
typedef void (APIENTRYP glVertexAttrib2f_t) (GLuint index, GLfloat x, GLfloat y);
static glVertexAttrib2f_t glVertexAttrib2f;
typedef void (APIENTRYP glVertexAttrib3f_t) (GLuint index, GLfloat x, GLfloat y, GLfloat z);
static glVertexAttrib3f_t glVertexAttrib3f;
typedef void (APIENTRYP glVertexAttrib4f_t) (GLuint index, GLfloat x, GLfloat y, GLfloat z, GLfloat w);
static glVertexAttrib4f_t glVertexAttrib4f;
typedef void (APIENTRYP glVertexAttribPointer_t) (GLuint index, GLint size, GLenum type, GLboolean normalized, GLsizei stride, const void *pointer);
static glVertexAttribPointer_t glVertexAttribPointer;
typedef void (APIENTRYP glEnableVertexAttribArray_t) (GLuint index);
static glEnableVertexAttribArray_t glEnableVertexAttribArray;
typedef void (APIENTRYP glDisableVertexAttribArray_t) (GLuint index);
static glDisableVertexAttribArray_t glDisableVertexAttribArray;
typedef void (APIENTRYP glGetActiveAttrib_t) (GLuint program, GLuint index, GLsizei bufSize, GLsizei *length, GLint *size, GLenum *type, GLchar *name);
static glGetActiveAttrib_t glGetActiveAttrib;
typedef GLint (APIENTRYP glGetAttribLocation_t) (GLuint program, const GLchar *name);
static glGetAttribLocation_t glGetAttribLocation;
// GL_VERSION_3_0
typedef void (APIENTRYP glGenVertexArrays_t) (GLsizei n, GLuint *arrays);
static glGenVertexArrays_t glGenVertexArrays;
typedef void (APIENTRYP glBindVertexArray_t) (GLuint array);
static glBindVertexArray_t glBindVertexArray;
typedef void (APIENTRYP glDeleteVertexArrays_t) (GLsizei n, const GLuint *arrays);
static glDeleteVertexArrays_t glDeleteVertexArrays;
typedef void (APIENTRYP glGenFramebuffers_t) (GLsizei n, GLuint *framebuffers);
static glGenFramebuffers_t glGenFramebuffers;
typedef void (APIENTRYP glFramebufferTexture2D_t) (GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level);
static glFramebufferTexture2D_t glFramebufferTexture2D;
typedef void (APIENTRYP glBindFramebuffer_t) (GLenum target, GLuint framebuffer);
static glBindFramebuffer_t glBindFramebuffer;
typedef void (APIENTRYP glDeleteFramebuffers_t) (GLsizei n, const GLuint *framebuffers);
static glDeleteFramebuffers_t glDeleteFramebuffers;
typedef void (APIENTRYP glClampColor_t) (GLenum target, GLenum clamp);
static glClampColor_t glClampColor;
// GL_ARB_color_buffer_float
typedef void (APIENTRYP glClampColorARB_t) (GLenum target, GLenum clamp);
static glClampColorARB_t glClampColorARB;
// GL_VERSION_3_1
typedef void (APIENTRYP glPrimitiveRestartIndex_t) (GLuint index);
static glPrimitiveRestartIndex_t glPrimitiveRestartIndex;
// GL_NV_primitive_restart
typedef void (APIENTRYP glPrimitiveRestartIndexNV_t) (GLuint index);
static glPrimitiveRestartIndexNV_t glPrimitiveRestartIndexNV;


// Load extension functions
typedef void (APIENTRYP GLFuncPtr)();

#if defined(_WIN32)

static GLFuncPtr glGetProcAddress(const char *name)
{
	return (GLFuncPtr) wglGetProcAddress((LPCSTR) name);
}

#elif defined(__APPLE__)

static GLFuncPtr glGetProcAddress(const char *name)
{
	static void *dlhandle = NULL;
	if (!dlhandle) {
		dlhandle = dlopen("/System/Library/Frameworks/OpenGL.framework/Versions/Current/OpenGL", RTLD_LAZY);
		if (!dlhandle)
			return NULL;
	}
	return (GLFuncPtr) dlsym(dlhandle, name);
}

#else

static GLFuncPtr glGetProcAddress(const char *name)
{
	return (GLFuncPtr) glXGetProcAddressARB((const GLubyte *) name);
}

#endif


// Load up all of the extension functions we use in this file
bool GLManager::glGetProcs()
{
	static bool already_tried = false;
	if (already_tried)
		return true;
	already_tried = true;

#define GETPROC(name, is_required) do { \
	if (!(name = (name ## _t) glGetProcAddress(#name)) && is_required) { \
		eprintf("GLManager::glGetProcs : couldn't get %s\n", #name); \
		return false; \
	} } while (0)

#ifdef _WIN32
	// GL_VERSION_1_3
	GETPROC(glActiveTexture, true);
#endif

	// GL_VERSION_1_5
	GETPROC(glGenBuffers, true);
	GETPROC(glBufferData, true);
	GETPROC(glBindBuffer, true);
	GETPROC(glDeleteBuffers, true);

	// GL_VERSION_2_0
	GETPROC(glCreateShader, true);
	GETPROC(glShaderSource, true);
	GETPROC(glCompileShader, true);
	GETPROC(glDeleteShader, true);
	GETPROC(glGetShaderiv, true);
	GETPROC(glGetShaderInfoLog, true);
	GETPROC(glCreateProgram, true);
	GETPROC(glAttachShader, true);
	GETPROC(glLinkProgram, true);
	GETPROC(glValidateProgram, true);
	GETPROC(glUseProgram, true);
	GETPROC(glDeleteProgram, true);
	GETPROC(glGetProgramiv, true);
	GETPROC(glGetProgramInfoLog, true);
	GETPROC(glUniform1f, true);
	GETPROC(glUniform2f, true);
	GETPROC(glUniform3f, true);
	GETPROC(glUniform4f, true);
	GETPROC(glUniform1i, true);
	GETPROC(glUniform2i, true);
	GETPROC(glUniform3i, true);
	GETPROC(glUniform4i, true);
	GETPROC(glUniformMatrix2fv, true);
	GETPROC(glUniformMatrix3fv, true);
	GETPROC(glUniformMatrix4fv, true);
	GETPROC(glGetActiveUniform, true);
	GETPROC(glGetUniformLocation, true);
	GETPROC(glVertexAttrib1f, true);
	GETPROC(glVertexAttrib2f, true);
	GETPROC(glVertexAttrib3f, true);
	GETPROC(glVertexAttrib4f, true);
	GETPROC(glVertexAttribPointer, true);
	GETPROC(glEnableVertexAttribArray, true);
	GETPROC(glDisableVertexAttribArray, true);
	GETPROC(glGetActiveAttrib, true);
	GETPROC(glGetAttribLocation, true);

	// GL_VERSION_3_0
	GETPROC(glGenVertexArrays, false);
	GETPROC(glBindVertexArray, false);
	GETPROC(glDeleteVertexArrays, false);
	GETPROC(glGenFramebuffers, false);
	GETPROC(glFramebufferTexture2D, false);
	GETPROC(glBindFramebuffer, false);
	GETPROC(glDeleteFramebuffers, false);
	GETPROC(glClampColor, false);

	// GL_ARB_color_buffer_float
	GETPROC(glClampColorARB, false);

	// GL_VERSION_3_1
	GETPROC(glPrimitiveRestartIndex, false);

	// GL_NV_primitive_restart
	GETPROC(glPrimitiveRestartIndexNV, false);

	return true;
#undef GETPROC
}


// Information about a #define to be passed to a shader.
class GLManager::ShaderDefine : public pair<string, string>
{
public:
	ShaderDefine(const string &defname, const string &value) :
		pair<string, string>(defname, value)
		{}
};


// Information about uniforms.
struct GLManager::UniformInfo
{
	string name;
	GLenum type;
	GLint location;
	UniformInfo() : type(GL_FLOAT), location(-1)
		{}
};


// Information about attributes.  We store this to emulate VAOs, plus
// legacy compatibility for old-style vertex buffer, etc.
struct GLManager::AttributeInfo
{
	string name;
	unsigned buf;
	const float *p;
	int floats_per_attribute;
	size_t stride, offset;
	GLint location;
	AttributeInfo() : buf(0), p(0),
		floats_per_attribute(0), stride(0), offset(0), location(-1)
		{}
};


// Information about a shader program and its components
class GLManager::ShaderInfo
{
public:
	string name;
	unsigned vshader, fshader, program;
	vector<bool> have_uniform, have_attribute;
	bool have_vertex, have_normal, have_color3, have_color4, have_texcoord;
	vector<UniformInfo> uniform_info;
	vector<AttributeInfo> attribute_info;
	AttributeInfo vertex_info;
	AttributeInfo normal_info;
	AttributeInfo color_info;
	AttributeInfo texcoord_info;
	vector<unsigned> texunit_binding;
	bool checked;

	ShaderInfo(const string &name_, unsigned vshader_, unsigned fshader_,
		unsigned program_, int nuniforms, int nattributes) :
		name(name_),
		vshader(vshader_), fshader(fshader_), program(program_),
		have_vertex(false), have_normal(false),
		have_color3(false), have_color4(false),
		have_texcoord(false), checked(false)
	{
		have_uniform.resize(nuniforms, false);
		uniform_info.resize(nuniforms);
		have_attribute.resize(nattributes, false);
		attribute_info.resize(nattributes);
		vertex_info.name = "gl_Vertex";
		vertex_info.floats_per_attribute = 3;
		normal_info.name = "gl_Normal";
		normal_info.floats_per_attribute = 3;
		color_info.name = "gl_Color";
		color_info.floats_per_attribute = 3;
		texcoord_info.name = "gl_MultiTexCoord0";
		texcoord_info.floats_per_attribute = 2;
	}
};


// Information about a VBO, IBO, texture, etc. - ID and where the data came from
class GLManager::BufferInfo
{
public:
	unsigned ind;
	const void *p;
	BufferInfo(unsigned ind_, const void *p_) : ind(ind_), p(p_)
		{}
};


// Set up state
GLManager::GLManager() : current_shader(0),
	current_buffer(0), current_ibuffer(0), current_texture(0),
	flipping_textures(false),
	checked_restart(false), have_restart_fixed(false),
	have_restart(false), have_restart_nv(false),
	checked_clamp(false), checked_vao(false), vao(0), eprintf_hook(NULL)
{
	// Compatibility: create default shader, buffer, and texture objects
	shaders.push_back(new ShaderInfo("No shader", 0, 0, 0, 0, 0));
	current_shader = shaders[0];
	shaders[0]->texunit_binding.push_back(0);
	buffers.push_back(new BufferInfo(0, 0));
	current_buffer = buffers.back();
	ibuffers.push_back(new BufferInfo(0, 0));
	current_ibuffer = ibuffers.back();
	textures.push_back(new BufferInfo(0, 0));
	current_texture = textures.back();
}


// Destructor - delete default objects
GLManager::~GLManager()
{
	clear();

	while (!shaders.empty()) {
		delete shaders.back();
		shaders.pop_back();
	}
	while (!buffers.empty()) {
		delete buffers.back();
		buffers.pop_back();
	}
	while (!ibuffers.empty()) {
		delete ibuffers.back();
		ibuffers.pop_back();
	}
	while (!textures.empty()) {
		delete textures.back();
		textures.pop_back();
	}

	// Delete VAO if we used one
	if (vao) {
		glBindVertexArray(0);
		glDeleteVertexArrays(1, &vao);
	}
}


// #define used by the next make_shader
void GLManager::define(const char *defname, const char *value)
{
	string defname_string(defname), value_string(value);

	// See if this define already exists.  If so, delete it.
	vector<ShaderDefine *>::iterator it = defines.begin();
	for ( ; it != defines.end(); it++) {
		if ((*it)->first == defname_string) {
			defines.erase(it);
			break;
		}
	}

	defines.push_back(new ShaderDefine(defname_string, value_string));
}


// Clear all defines.  Note that make_shader() does *not* do this.
void GLManager::clear_defines()
{
	while (!defines.empty()) {
		delete defines.back();
		defines.pop_back();
	}
}


// Create a shader program using the given vertex and fragment shaders
unsigned GLManager::make_shader(const char *name,
	const char *vert_text, const char *frag_text)
{
	// Sanity checking
	if (!name || !*name) {
		eprintf("GLManager::make_shader : empty name\n");
		return 0;
	}
	if (!vert_text || !*vert_text) {
		eprintf("GLManager::make_shader : empty vertex shader\n");
		return 0;
	}
	if (!frag_text || !*frag_text) {
		eprintf("GLManager::make_shader : empty fragment shader\n");
		return 0;
	}
	if (!glGetProcs())
		return 0;

	// Build strings containing defines plus shader code
	string defines_string;
	for (size_t i = 0; i < defines.size(); i++) {
		defines_string += "#define " + defines[i]->first + " "
			+ defines[i]->second + "\n";
	}
	defines_string += "#line 0\n";
	string vert_string = defines_string + vert_text;
	string frag_string = defines_string + frag_text;

	// Compile and link
	unsigned vshader = glCreateShader(GL_VERTEX_SHADER);
	const char *vsource = vert_string.c_str();
	glShaderSource(vshader, 1, &vsource, NULL);
	glCompileShader(vshader);
	if (!shader_ok(vshader, "compiling vertex shader", name)) {
		glDeleteShader(vshader);
		return 0;
	}

	unsigned fshader = glCreateShader(GL_FRAGMENT_SHADER);
	const char *fsource = frag_string.c_str();
	glShaderSource(fshader, 1, &fsource, NULL);
	glCompileShader(fshader);
	if (!shader_ok(fshader, "compiling fragment shader", name)) {
		glDeleteShader(fshader);
		glDeleteShader(vshader);
		return 0;
	}

	unsigned program = glCreateProgram();
	glAttachShader(program, vshader);
	glAttachShader(program, fshader);
	glLinkProgram(program);
	if (!program_ok(program, "linking shader program", name)) {
		glDeleteProgram(program);
		glDeleteShader(fshader);
		glDeleteShader(vshader);
		return 0;
	}

	// Figure out number of uniforms and attributes
	int nuniforms = 0, nattributes = 0;
	glGetProgramiv(program, GL_ACTIVE_UNIFORMS, &nuniforms);
	glGetProgramiv(program, GL_ACTIVE_ATTRIBUTES, &nattributes);

	// Create ShaderInfo node
	shaders.push_back(new ShaderInfo(name, vshader, fshader, program,
		nuniforms, nattributes));

	// Fill in UniformInfo and AttributeInfo nodes
	int len = 1024;
	// Buggy on some drivers...
	//glGetProgramiv(program, GL_ACTIVE_UNIFORM_MAX_LENGTH, &len);
	vector<char> uname(len);
	for (int i = 0; i < nuniforms; i++) {
		GLint size;
		GLenum type;
		glGetActiveUniform(program, i, len, NULL,
			&size, &type, &uname[0]);
		UniformInfo &ui = shaders.back()->uniform_info[i];
		ui.name = string(&uname[0]);
		ui.type = type;
		if (type == GL_SAMPLER_2D)
			shaders.back()->texunit_binding.push_back(0);
		ui.location = glGetUniformLocation(program, &uname[0]);
	}

	// Buggy on some drivers...
	//glGetProgramiv(program, GL_ACTIVE_ATTRIBUTE_MAX_LENGTH, &len);
	vector<char> aname(len);
	for (int i = 0; i < nattributes; i++) {
		GLint size;
		GLenum type;
		glGetActiveAttrib(program, i, len, NULL,
			&size, &type, &aname[0]);
		AttributeInfo &ai = shaders.back()->attribute_info[i];
		ai.name = string(&aname[0]);
		if (type == GL_FLOAT)
			ai.floats_per_attribute = size;
		else if (type == GL_FLOAT_VEC2)
			ai.floats_per_attribute = 2 * size;
		else if (type == GL_FLOAT_VEC3)
			ai.floats_per_attribute = 3 * size;
		else if (type == GL_FLOAT_VEC4)
			ai.floats_per_attribute = 4 * size;
		ai.location = glGetAttribLocation(program, &aname[0]);
	}

	// Make this shader active
//	check_vao();
	use_shader(program);
	return program;
}


// Create a shader program by reading shaders from files
unsigned GLManager::make_shader_from_files(const char *name,
	const char *vert_filename, const char *frag_filename)
{
	ifstream vert_in(vert_filename, ios::in | ios::binary);
	if (!vert_in) {
		eprintf("GLManager::make_shader_from_files : couldn't open %s\n", vert_filename);
		return 0;
	}
	ifstream frag_in(frag_filename, ios::in | ios::binary);
	if (!frag_in) {
		eprintf("GLManager::make_shader_from_files : couldn't open %s\n", frag_filename);
		return 0;
	}

	string vert_string, frag_string;
	vert_in.seekg(0, ios::end);
	vert_string.resize(vert_in.tellg());
	vert_in.seekg(0, ios::beg);
	vert_in.read(&vert_string[0], vert_string.size());
	vert_in.close();
	frag_in.seekg(0, ios::end);
	frag_string.resize(frag_in.tellg());
	frag_in.seekg(0, ios::beg);
	frag_in.read(&frag_string[0], frag_string.size());
	frag_in.close();

	return make_shader(name, vert_string.c_str(), frag_string.c_str());
}


// Make shader active, for rendering and uniform/attribute definitions
bool GLManager::use_shader(unsigned ind)
{
	if (current_shader->program == ind)
		return true;

	// Disable old attribute arrays
	if (current_shader->have_vertex)
		glDisableClientState(GL_VERTEX_ARRAY);
	if (current_shader->have_normal)
		glDisableClientState(GL_NORMAL_ARRAY);
	if (current_shader->have_color3 || current_shader->have_color4)
		glDisableClientState(GL_COLOR_ARRAY);
	if (current_shader->have_texcoord)
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	for (size_t i = 0; i < current_shader->attribute_info.size(); i++)
		if (current_shader->have_attribute[i])
			glDisableVertexAttribArray(current_shader->attribute_info[i].location);

	// Find new shader
	for (size_t i = 0; i < shaders.size(); i++) {
		if (shaders[i]->program == ind) {
			current_shader = shaders[i];
			break;
		}
	}
	if (current_shader->program != ind) {
		use_shader();
		eprintf("GLManager::use_shader : can't find shader program %u\n", ind);
		return false;
	}
	glUseProgram(current_shader->program);

	// For compatibility, simulate the effect of VAOs
	if (current_shader->have_vertex) {
		AttributeInfo &ai = current_shader->vertex_info;
		if (ai.buf)
			vertexarray3fv(ai.buf, ai.stride, ai.offset);
		else if (ai.p)
			vertexarray3fv(ai.p, ai.stride, ai.offset);
	}
	if (current_shader->have_normal) {
		AttributeInfo &ai = current_shader->normal_info;
		if (ai.buf)
			normalarray3fv(ai.buf, ai.stride, ai.offset);
		else if (ai.p)
			normalarray3fv(ai.p, ai.stride, ai.offset);
	}
	if (current_shader->have_color3) {
		AttributeInfo &ai = current_shader->color_info;
		if (ai.buf)
			colorarray3fv(ai.buf, ai.stride, ai.offset);
		else if (ai.p)
			colorarray3fv(ai.p, ai.stride, ai.offset);
	}
	if (current_shader->have_color4) {
		AttributeInfo &ai = current_shader->color_info;
		if (ai.buf)
			colorarray4fv(ai.buf, ai.stride, ai.offset);
		else if (ai.p)
			colorarray4fv(ai.p, ai.stride, ai.offset);
	}
	if (current_shader->have_texcoord) {
		AttributeInfo &ai = current_shader->texcoord_info;
		if (ai.buf)
			texcoordarray2fv(ai.buf, ai.stride, ai.offset);
		else if (ai.p)
			texcoordarray2fv(ai.p, ai.stride, ai.offset);
	}
	for (size_t i = 0; i < current_shader->attribute_info.size(); i++) {
		if (!current_shader->have_attribute[i])
			continue;
		AttributeInfo &ai = current_shader->attribute_info[i];
		if (ai.buf)
			attributearray(ai.location, ai.floats_per_attribute, ai.buf, ai.stride, ai.offset);
		else if (ai.p)
			attributearray(ai.location, ai.floats_per_attribute, ai.p, ai.stride, ai.offset);
	}

	// Bind the textures we're using
	for (size_t i = 0; i < current_shader->texunit_binding.size(); i++) {
		unsigned texind = current_shader->texunit_binding[i];
		glActiveTexture(GL_TEXTURE0 + i);
		glBindTexture(GL_TEXTURE_2D, texind);
		if (ind == 0) {
			if (texind)
				glEnable(GL_TEXTURE_2D);
			else
				glDisable(GL_TEXTURE_2D);
		}
	}
	glActiveTexture(GL_TEXTURE0);
	return true;
}


// Delete a shader
void GLManager::clear_shader(unsigned ind)
{
	if (current_shader->program == ind)
		use_shader();

	vector<ShaderInfo *>::iterator it = shaders.begin();
	for (it++; it != shaders.end(); it++) {
		if ((*it)->program == ind)
			break;
	}
	if (it == shaders.end()) {
		eprintf("GLManager::clear_shader : can't find shader program %u\n", ind);
		return;
	}

	if ((*it)->program)
		glDeleteProgram((*it)->program);
	if ((*it)->fshader)
		glDeleteShader((*it)->fshader);
	if ((*it)->vshader)
		glDeleteShader((*it)->vshader);
	delete (*it);
	shaders.erase(it);
}


// Delete all shaders
void GLManager::clear_shaders()
{
	use_shader();
	while (shaders.size() > 1)
		clear_shader(shaders.back()->program);
}


// Set a uniform for the current shader
bool GLManager::uniform_boilerplate(int ind, unsigned type)
{
	if (!current_shader->program) {
		eprintf("GLManager::uniform : no active shader program\n");
		return false;
	}
	size_t which = 0;
	for (which = 0; which < current_shader->uniform_info.size(); which++)
		if (current_shader->uniform_info[which].location == ind)
			break;
	if (which == current_shader->uniform_info.size()) {
		eprintf("GLManager::uniform : invalid uniform\n");
		return false;
	}
	if (current_shader->uniform_info[which].type != (GLenum) type) {
		eprintf("GLManager::uniform : type mismatch specifying uniform %s\n",
			current_shader->uniform_info[which].name.c_str());
		return false;
	}
	current_shader->have_uniform[which] = true;
	return true;
}

void GLManager::uniform1i(int ind, int x)
	{ if (uniform_boilerplate(ind, GL_INT)) glUniform1i(ind, x); }
void GLManager::uniform2i(int ind, int x, int y)
	{ if (uniform_boilerplate(ind, GL_INT_VEC2)) glUniform2i(ind, x, y); }
void GLManager::uniform3i(int ind, int x, int y, int z)
	{ if (uniform_boilerplate(ind, GL_INT_VEC3)) glUniform3i(ind, x, y, z); }
void GLManager::uniform4i(int ind, int x, int y, int z, int w)
	{ if (uniform_boilerplate(ind, GL_INT_VEC4)) glUniform4i(ind, x, y, z, w); }
void GLManager::uniform1f(int ind, float x)
	{ if (uniform_boilerplate(ind, GL_FLOAT)) glUniform1f(ind, x); }
void GLManager::uniform2f(int ind, float x, float y)
	{ if (uniform_boilerplate(ind, GL_FLOAT_VEC2)) glUniform2f(ind, x, y); }
void GLManager::uniform3f(int ind, float x, float y, float z)
	{ if (uniform_boilerplate(ind, GL_FLOAT_VEC3)) glUniform3f(ind, x, y, z); }
void GLManager::uniform4f(int ind, float x, float y, float z, float w)
	{ if (uniform_boilerplate(ind, GL_FLOAT_VEC4)) glUniform4f(ind, x, y, z, w); }
void GLManager::uniformMatrix2fv(int ind, const float *p)
	{ if (uniform_boilerplate(ind, GL_FLOAT_MAT2)) glUniformMatrix2fv(ind, 1, false, p); }
void GLManager::uniformMatrix3fv(int ind, const float *p)
	{ if (uniform_boilerplate(ind, GL_FLOAT_MAT3)) glUniformMatrix3fv(ind, 1, false, p); }
void GLManager::uniformMatrix4fv(int ind, const float *p)
	{ if (uniform_boilerplate(ind, GL_FLOAT_MAT4)) glUniformMatrix4fv(ind, 1, false, p); }

void GLManager::uniformTexture(int ind, unsigned texind)
{
	if (!uniform_boilerplate(ind, GL_SAMPLER_2D))
		return;
	unsigned texunit = 0;
	for (int i = 0; i < ind; i++) {
		if (current_shader->uniform_info[i].type == GL_SAMPLER_2D)
			texunit++;
	}

	glActiveTexture(GL_TEXTURE0 + texunit);
	use_texture(texind);
	glUniform1i(ind, int(texunit));
}


// "Clear" uniforms so that they have to be set again
void GLManager::clear_uniform(int ind)
{
	if (ind < 0 || ind >= (int) current_shader->have_uniform.size()) {
		eprintf("GLManager::clear_uniform : invalid uniform\n");
		return;
	}
	current_shader->have_uniform[ind] = false;
	current_shader->checked = false;
}

void GLManager::clear_uniforms()
{
	for (size_t i = 0; i < current_shader->have_uniform.size(); i++)
		current_shader->have_uniform[i] = false;
	current_shader->checked = false;
}


// Set an attribute to a (client-side) array or vector
void GLManager::attributearray(int ind, int floats_per_attribute,
	const float *p, size_t stride, size_t offset)
{
	if (!current_shader->program) {
		eprintf("GLManager::attribute : no active shader program\n");
		return;
	}
	size_t which = 0;
	for (which = 0; which < current_shader->attribute_info.size(); which++)
		if (current_shader->attribute_info[which].location == ind)
			break;
	if (which == current_shader->attribute_info.size()) {
		eprintf("GLManager::attribute : invalid attribute\n");
		return;
	}
	if (current_shader->attribute_info[which].floats_per_attribute !=
		floats_per_attribute) {
		eprintf("GLManager::attribute : type mismatch specifying attribute %s\n",
			current_shader->attribute_info[which].name.c_str());
		return;
	}
	use_buffer();
	current_shader->have_attribute[which] = true;
	AttributeInfo &ai = current_shader->attribute_info[which];
	ai.buf = 0;
	ai.p = p;
	ai.stride = stride;
	ai.offset = offset;
	glVertexAttribPointer(ind, floats_per_attribute, GL_FLOAT, GL_FALSE,
		stride, p + offset / sizeof(float));
	glEnableVertexAttribArray(ind);
}


// Set an attribute to a buffer
void GLManager::attributearray(int ind, int floats_per_attribute, unsigned buf,
	size_t stride, size_t offset)
{
	if (!current_shader->program) {
		eprintf("GLManager::attribute : no active shader program\n");
		return;
	}
	size_t which = 0;
	for (which = 0; which < current_shader->attribute_info.size(); which++)
		if (current_shader->attribute_info[which].location == ind)
			break;
	if (which == current_shader->attribute_info.size()) {
		eprintf("GLManager::attribute : invalid attribute\n");
		return;
	}
	AttributeInfo &ai = current_shader->attribute_info[which];
	if (ai.floats_per_attribute != floats_per_attribute) {
		eprintf("GLManager::attribute : type mismatch specifying attribute %s\n",
			current_shader->attribute_info[which].name.c_str());
		return;
	}
	if (!use_buffer(buf))
		return;
	current_shader->have_attribute[which] = true;
	ai.buf = buf;
	ai.p = 0;
	ai.stride = stride;
	ai.offset = offset;
	glVertexAttribPointer(ind, floats_per_attribute, GL_FLOAT, GL_FALSE,
		stride, GL_VBO_OFFSET(offset));
	glEnableVertexAttribArray(ind);
}


// Set an attribute to a constant
bool GLManager::attribute_boilerplate(int ind, int floats_per_attribute)
{
	if (!current_shader->program) {
		eprintf("GLManager::attribute : no active shader program\n");
		return false;
	}
	size_t which = 0;
	for (which = 0; which < current_shader->attribute_info.size(); which++)
		if (current_shader->attribute_info[which].location == ind)
			break;
	if (which == current_shader->attribute_info.size()) {
		eprintf("GLManager::attribute : invalid attribute\n");
		return false;
	}
	AttributeInfo &ai = current_shader->attribute_info[which];
	if (ai.floats_per_attribute != floats_per_attribute) {
		eprintf("GLManager::attribute : type mismatch specifying attribute %s\n",
			current_shader->attribute_info[which].name.c_str());
		return false;
	}
	current_shader->have_attribute[which] = true;
	ai.buf = 0;
	ai.p = 0;
	ai.stride = 0;
	ai.offset = 0;
	glDisableVertexAttribArray(ind);
	return true;
}

void GLManager::attribute1f(int ind, float x)
	{ if (attribute_boilerplate(ind, 1)) glVertexAttrib1f(ind, x); }
void GLManager::attribute2f(int ind, float x, float y)
	{ if (attribute_boilerplate(ind, 2)) glVertexAttrib2f(ind, x, y); }
void GLManager::attribute3f(int ind, float x, float y, float z)
	{ if (attribute_boilerplate(ind, 3)) glVertexAttrib3f(ind, x, y, z); }
void GLManager::attribute4f(int ind, float x, float y, float z, float w)
	{ if (attribute_boilerplate(ind, 4)) glVertexAttrib4f(ind, x, y, z, w); }


// Clear attribute bindings
void GLManager::clear_attribute(int ind)
{
	if (ind < 0 || ind >= (int) current_shader->have_attribute.size()) {
		eprintf("GLManager::clear_attribute : invalid attribute\n");
		return;
	}
	current_shader->have_attribute[ind] = false;
	glDisableVertexAttribArray(ind);
	current_shader->checked = false;
}

void GLManager::clear_attributes()
{
	clear_vertexarray();
	clear_normalarray();
	clear_colorarray();
	clear_texcoordarray();
	for (size_t i = 0; i < current_shader->have_attribute.size(); i++) {
		current_shader->have_attribute[i] = false;
		glDisableVertexAttribArray(i);
	}
	current_shader->checked = false;
}


// Old-style attributes
void GLManager::vertexarray3fv(const float *p, size_t stride, size_t offset)
{
	use_buffer();
	current_shader->have_vertex = true;
	AttributeInfo &ai = current_shader->vertex_info;
	ai.buf = 0;
	ai.p = p;
	ai.stride = stride;
	ai.offset = offset;
	glVertexPointer(3, GL_FLOAT, stride, p + offset / sizeof(float));
	glEnableClientState(GL_VERTEX_ARRAY);
}

void GLManager::normalarray3fv(const float *p, size_t stride, size_t offset)
{
	use_buffer();
	current_shader->have_normal = true;
	AttributeInfo &ai = current_shader->normal_info;
	ai.buf = 0;
	ai.p = p;
	ai.stride = stride;
	ai.offset = offset;
	glNormalPointer(GL_FLOAT, stride, p + offset / sizeof(float));
	glEnableClientState(GL_NORMAL_ARRAY);
}

void GLManager::colorarray3fv(const float *p, size_t stride, size_t offset)
{
	use_buffer();
	current_shader->have_color3 = true;
	current_shader->have_color4 = false;
	AttributeInfo &ai = current_shader->color_info;
	ai.buf = 0;
	ai.p = p;
	ai.stride = stride;
	ai.offset = offset;
	glColorPointer(3, GL_FLOAT, stride, p + offset / sizeof(float));
	glEnableClientState(GL_COLOR_ARRAY);
}

void GLManager::colorarray4fv(const float *p, size_t stride, size_t offset)
{
	use_buffer();
	current_shader->have_color3 = false;
	current_shader->have_color4 = true;
	AttributeInfo &ai = current_shader->color_info;
	ai.buf = 0;
	ai.p = p;
	ai.stride = stride;
	ai.offset = offset;
	glColorPointer(4, GL_FLOAT, stride, p + offset / sizeof(float));
	glEnableClientState(GL_COLOR_ARRAY);
}

void GLManager::texcoordarray2fv(const float *p, size_t stride, size_t offset)
{
	use_buffer();
	current_shader->have_texcoord = true;
	AttributeInfo &ai = current_shader->texcoord_info;
	ai.buf = 0;
	ai.p = p;
	ai.stride = stride;
	ai.offset = offset;
	glTexCoordPointer(2, GL_FLOAT, stride, p + offset / sizeof(float));
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
}

void GLManager::vertexarray3fv(unsigned buf, size_t stride, size_t offset)
{
	if (!use_buffer(buf))
		return;
	current_shader->have_vertex = true;
	AttributeInfo &ai = current_shader->vertex_info;
	ai.buf = buf;
	ai.p = 0;
	ai.stride = stride;
	ai.offset = offset;
	glVertexPointer(3, GL_FLOAT, stride, GL_VBO_OFFSET(offset));
	glEnableClientState(GL_VERTEX_ARRAY);
}

void GLManager::normalarray3fv(unsigned buf, size_t stride, size_t offset)
{
	if (!use_buffer(buf))
		return;
	current_shader->have_normal = true;
	AttributeInfo &ai = current_shader->normal_info;
	ai.buf = buf;
	ai.p = 0;
	ai.stride = stride;
	ai.offset = offset;
	glNormalPointer(GL_FLOAT, stride, GL_VBO_OFFSET(offset));
	glEnableClientState(GL_NORMAL_ARRAY);
}

void GLManager::colorarray3fv(unsigned buf, size_t stride, size_t offset)
{
	if (!use_buffer(buf))
		return;
	current_shader->have_color3 = true;
	current_shader->have_color4 = false;
	AttributeInfo &ai = current_shader->color_info;
	ai.buf = buf;
	ai.p = 0;
	ai.stride = stride;
	ai.offset = offset;
	glColorPointer(3, GL_FLOAT, stride, GL_VBO_OFFSET(offset));
	glEnableClientState(GL_COLOR_ARRAY);
}

void GLManager::colorarray4fv(unsigned buf, size_t stride, size_t offset)
{
	if (!use_buffer(buf))
		return;
	current_shader->have_color3 = false;
	current_shader->have_color4 = true;
	AttributeInfo &ai = current_shader->color_info;
	ai.buf = buf;
	ai.p = 0;
	ai.stride = stride;
	ai.offset = offset;
	glColorPointer(4, GL_FLOAT, stride, GL_VBO_OFFSET(offset));
	glEnableClientState(GL_COLOR_ARRAY);
}

void GLManager::texcoordarray2fv(unsigned buf, size_t stride, size_t offset)
{
	if (!use_buffer(buf))
		return;
	current_shader->have_texcoord = true;
	AttributeInfo &ai = current_shader->texcoord_info;
	ai.buf = buf;
	ai.p = 0;
	ai.stride = stride;
	ai.offset = offset;
	glTexCoordPointer(2, GL_FLOAT, stride, GL_VBO_OFFSET(offset));
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
}

void GLManager::clear_vertexarray()
{
	if (current_shader->have_vertex) {
		glDisableClientState(GL_VERTEX_ARRAY);
		current_shader->have_vertex = false;
	}
}

void GLManager::clear_normalarray()
{
	if (current_shader->have_normal) {
		glDisableClientState(GL_NORMAL_ARRAY);
		current_shader->have_normal = false;
	}
}

void GLManager::clear_colorarray()
{
	if (current_shader->have_color3 || current_shader->have_color4) {
		glDisableClientState(GL_COLOR_ARRAY);
		current_shader->have_color3 = false;
		current_shader->have_color4 = false;
	}
}

void GLManager::clear_texcoordarray()
{
	if (current_shader->have_texcoord) {
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		current_shader->have_texcoord = false;
	}
}


// Old-style constant attributes
void GLManager::normal3f(float nx, float ny, float nz)
{
	current_shader->have_normal = true;
	AttributeInfo &ai = current_shader->normal_info;
	ai.buf = 0;
	ai.p = 0;
	glDisableClientState(GL_NORMAL_ARRAY);
	glNormal3f(nx, ny, nz);
}

void GLManager::color3f(float r, float g, float b)
{
	current_shader->have_color3 = true;
	current_shader->have_color4 = false;
	AttributeInfo &ai = current_shader->color_info;
	ai.buf = 0;
	ai.p = 0;
	glDisableClientState(GL_COLOR_ARRAY);
	glColor3f(r, g, b);
}

void GLManager::color4f(float r, float g, float b, float a)
{
	current_shader->have_color3 = false;
	current_shader->have_color4 = true;
	AttributeInfo &ai = current_shader->color_info;
	ai.buf = 0;
	ai.p = 0;
	glDisableClientState(GL_COLOR_ARRAY);
	glColor4f(r, g, b, a);
}

void GLManager::texcoord2f(float u, float v)
{
	current_shader->have_texcoord = true;
	AttributeInfo &ai = current_shader->texcoord_info;
	ai.buf = 0;
	ai.p = 0;
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glTexCoord2f(u, v);
}


// Create a buffer and copy in some data
unsigned GLManager::make_buffer(const float *p, size_t nfloat)
{
	if (!glGetProcs())
		return 0;
	unsigned buf;
	glGenBuffers(1, &buf);
	buffers.push_back(new BufferInfo(buf, p));
	update_buffer(buf, p, nfloat);
	return buf;
}

void GLManager::update_buffer(unsigned buf, const float *p, size_t nfloat)
{
	if (!use_buffer(buf))
		return;
	glBufferData(GL_ARRAY_BUFFER, nfloat * sizeof(float), p, GL_STATIC_DRAW);
	for (size_t i = 1; i < buffers.size(); i++) {
		if (buffers[i]->ind == buf) {
			buffers[i]->p = p;
			return;
		}
	}
}

bool GLManager::use_buffer(unsigned buf)
{
	if (current_buffer->ind == buf)
		return true;
	if (!glGetProcs())
		return false;
	for (size_t i = 0; i < buffers.size(); i++) {
		if (buffers[i]->ind == buf) {
			current_buffer = buffers[i];
			break;
		}
	}
	if (current_buffer->ind != buf) {
		use_buffer();
		eprintf("GLManager::use_buffer : can't find buffer %u\n", buf);
		return false;
	}
	glBindBuffer(GL_ARRAY_BUFFER, buf);
	return true;
}

void GLManager::clear_buffer(unsigned buf)
{
	if (current_buffer->ind == buf)
		use_buffer();

	vector<BufferInfo *>::iterator it = buffers.begin();
	for (it++; it != buffers.end(); it++) {
		if ((*it)->ind == buf)
			break;
	}
	if (it == buffers.end()) {
		eprintf("GLManager::clear_buffer : can't find buffer %u\n", buf);
		return;
	}

	if ((*it)->ind)
		glDeleteBuffers(1, &((*it)->ind));
	delete (*it);
	buffers.erase(it);
}

void GLManager::clear_buffers()
{
	use_buffer();
	while (buffers.size() > 1)
		clear_buffer(buffers.back()->ind);
}

unsigned GLManager::make_ibuffer(const unsigned *p, size_t ninds)
{
	if (!glGetProcs())
		return 0;
	unsigned buf;
	glGenBuffers(1, &buf);
	ibuffers.push_back(new BufferInfo(buf, p));
	update_ibuffer(buf, p, ninds);
	return buf;
}

void GLManager::update_ibuffer(unsigned buf, const unsigned *p, size_t ninds)
{
	if (!use_ibuffer(buf))
		return;
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, ninds * sizeof(unsigned), p, GL_STATIC_DRAW);
	for (size_t i = 1; i < ibuffers.size(); i++) {
		if (ibuffers[i]->ind == buf) {
			ibuffers[i]->p = p;
			return;
		}
	}
}

bool GLManager::use_ibuffer(unsigned buf)
{
	if (current_ibuffer->ind == buf)
		return true;
	if (!glGetProcs())
		return false;
	for (size_t i = 0; i < ibuffers.size(); i++) {
		if (ibuffers[i]->ind == buf) {
			current_ibuffer = ibuffers[i];
			break;
		}
	}
	if (current_ibuffer->ind != buf) {
		use_ibuffer();
		eprintf("GLManager::use_ibuffer : can't find buffer %u\n", buf);
		return false;
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf);
	return true;
}

void GLManager::clear_ibuffer(unsigned buf)
{
	if (current_ibuffer->ind == buf)
		use_ibuffer();

	vector<BufferInfo *>::iterator it = ibuffers.begin();
	for (it++; it != ibuffers.end(); it++) {
		if ((*it)->ind == buf)
			break;
	}
	if (it == ibuffers.end()) {
		eprintf("GLManager::clear_ibuffer : can't find buffer %u\n", buf);
		return;
	}

	if ((*it)->ind)
		glDeleteBuffers(1, &((*it)->ind));
	delete (*it);
	ibuffers.erase(it);
}

void GLManager::clear_ibuffers()
{
	use_ibuffer();
	while (ibuffers.size() > 1)
		clear_ibuffer(ibuffers.back()->ind);
}


// Length of a pixel in bytes
int GLManager::pixel_in_bytes(unsigned format, unsigned type)
{
	int components = 1;
	switch (format) {
		case GL_RG:
		case GL_DEPTH_STENCIL:
		case GL_LUMINANCE_ALPHA:
			components = 2;
			break;
		case GL_RGB:
		case GL_BGR:
			components = 3;
			break;
		case GL_RGBA:
		case GL_BGRA:
			components = 4;
			break;
	}

	int bytes = components;
	switch (type) {
		case GL_UNSIGNED_SHORT:
		case GL_SHORT:
		case GL_HALF_FLOAT:
			bytes = 2 * components;
			break;
		case GL_UNSIGNED_INT:
		case GL_INT:
		case GL_FLOAT:
			bytes = 4 * components;
			break;
		case GL_UNSIGNED_BYTE_3_3_2:
		case GL_UNSIGNED_BYTE_2_3_3_REV:
			bytes = 1;
			break;
		case GL_UNSIGNED_SHORT_5_6_5:
		case GL_UNSIGNED_SHORT_5_6_5_REV:
		case GL_UNSIGNED_SHORT_4_4_4_4:
		case GL_UNSIGNED_SHORT_4_4_4_4_REV:
		case GL_UNSIGNED_SHORT_5_5_5_1:
		case GL_UNSIGNED_SHORT_1_5_5_5_REV:
			bytes = 2;
			break;
		case GL_UNSIGNED_INT_8_8_8_8:
		case GL_UNSIGNED_INT_8_8_8_8_REV:
		case GL_UNSIGNED_INT_10_10_10_2:
		case GL_UNSIGNED_INT_2_10_10_10_REV:
			bytes = 4;
			break;
	}
	return bytes;
}


// Should future texture loads flip the image in Y?
void GLManager::flip_textures(bool flip)
{
	flipping_textures = flip;
}


// Textures
unsigned GLManager::make_texture(int w, int h, const void *p,
	unsigned format, unsigned type, unsigned internalformat)
{
	if (!glGetProcs())
		return 0;
	unsigned texind;
	glGenTextures(1, &texind);
	textures.push_back(new BufferInfo(texind, p));
	update_texture(texind, w, h, p, format, type, internalformat);
	return texind;
}

void GLManager::update_texture(unsigned texind, int w, int h, const void *p,
	unsigned format, unsigned type, unsigned internalformat)
{
	if (!use_texture(texind))
		return;
	if (!internalformat)
		internalformat = format;

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	if (flipping_textures) {
		int rowlen = w * pixel_in_bytes(format, type);
		unsigned char *buf = new unsigned char[rowlen * h];
		for (int i = 0; i < h; i++) {
			memcpy(buf + rowlen * (h - 1 - i),
			       (const unsigned char *) p + rowlen * i,
			       rowlen);
		}
		glTexImage2D(GL_TEXTURE_2D, 0, internalformat, w, h, 0, format, type, buf);
		delete [] buf;
	} else {
		glTexImage2D(GL_TEXTURE_2D, 0, internalformat, w, h, 0, format, type, p);
	}

	glBindTexture(GL_TEXTURE_2D, texind);  // Work around some driver bugs...
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	if (internalformat == GL_R16F || internalformat == GL_R32F ||
	    internalformat == GL_RG16F || internalformat == GL_RG32F ||
	    internalformat == GL_RGB16F || internalformat == GL_RGB32F ||
	    internalformat == GL_RGBA16F || internalformat == GL_RGBA32F) {
		check_clamp();
	}

	for (size_t i = 1; i < textures.size(); i++) {
		if (textures[i]->ind == texind) {
			textures[i]->p = p;
			return;
		}
	}
}

bool GLManager::use_texture(unsigned texind)
{
	if (!glGetProcs())
		return false;
	unsigned texunit = 0;
	glGetIntegerv(GL_ACTIVE_TEXTURE, (GLint *) &texunit);
	texunit -= GL_TEXTURE0;

	if (current_texture->ind == texind &&
	    texunit < current_shader->texunit_binding.size() &&
	    current_shader->texunit_binding[texunit] == texind)
		return true;

	for (size_t i = 0; i < textures.size(); i++) {
		if (textures[i]->ind == texind) {
			current_texture = textures[i];
			break;
		}
	}
	if (current_texture->ind != texind) {
		use_texture();
		eprintf("GLManager::use_texture : can't find texture %u\n", texind);
		return false;
	}

	glBindTexture(GL_TEXTURE_2D, texind);
	if (texunit < current_shader->texunit_binding.size())
		current_shader->texunit_binding[texunit] = texind;
	if (current_shader->program == 0) {
		if (texind)
			glEnable(GL_TEXTURE_2D);
		else
			glDisable(GL_TEXTURE_2D);
	}
	return true;
}

void GLManager::clear_texture(unsigned texind)
{
	if (current_texture->ind == texind)
		use_texture();

	vector<BufferInfo *>::iterator it = textures.begin();
	for (it++; it != textures.end(); it++) {
		if ((*it)->ind == texind)
			break;
	}
	if (it == textures.end()) {
		eprintf("GLManager::clear_texture : can't find texture %u\n", texind);
		return;
	}

	if ((*it)->ind)
		glDeleteTextures(1, &((*it)->ind));
	delete (*it);
	textures.erase(it);
}

void GLManager::clear_textures()
{
	use_texture();
	while (textures.size() > 1)
		clear_texture(textures.back()->ind);
}


// Render to texture
void GLManager::rendertarget(unsigned texind, unsigned depthbufferind)
{
	if (!glGetProcs() || glversion() < 3.0f)
		return;

	unsigned fb;
	glGetIntegerv(GL_FRAMEBUFFER_BINDING, (GLint *) &fb);
	if (texind == 0) {
		if (fb != 0) {
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			glDeleteFramebuffers(1, &fb);
			glDrawBuffer(GL_BACK);
		}
		return;
	}

	if (fb == 0) {
		glGenFramebuffers(1, &fb);
		glBindFramebuffer(GL_FRAMEBUFFER, fb);
	}

	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
		GL_TEXTURE_2D, texind, 0);
	if (depthbufferind)
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
			GL_TEXTURE_2D, depthbufferind, 0);
	glDrawBuffer(GL_COLOR_ATTACHMENT0);
}


// Is shader program ready to go, and have all required
// attributes and uniforms been specified?
bool GLManager::shader_ready()
{
	int program = current_shader->program;
	if (!program || current_shader->checked)
		return true;
	current_shader->checked = true;

	// Validate
	glValidateProgram(program);
	GLint status, len;
	glGetProgramiv(program, GL_VALIDATE_STATUS, &status);
	if (status != GL_TRUE) {
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
		vector<char> log(len);
		glGetProgramInfoLog(program, len, NULL, &log[0]);

		eprintf("GLManager::shader_ready : %s failed to validate\n%s\n",
			current_shader->name.c_str(), &log[0]);
		return false;
	}

	bool ready = true;

	// Check uniforms
	for (size_t i = 0; i < current_shader->have_uniform.size(); i++) {
		if (current_shader->uniform_info[i].location < 0 ||
		    current_shader->have_uniform[i])
			continue;
		eprintf("GLManager::shader_ready : shader %s failed to set uniform %s\n",
			current_shader->name.c_str(),
			current_shader->uniform_info[i].name.c_str());
		ready = false;
	}

	// Check attributes
	for (size_t i = 0; i < current_shader->have_attribute.size(); i++) {
		if (current_shader->attribute_info[i].location < 0 ||
		    current_shader->have_attribute[i])
			continue;
		const char *aname = current_shader->attribute_info[i].name.c_str();
		if (strcmp(aname, current_shader->vertex_info.name.c_str()) == 0 && current_shader->have_vertex)
			continue;
		if (strcmp(aname, current_shader->normal_info.name.c_str()) == 0 && current_shader->have_normal)
			continue;
		if (strcmp(aname, current_shader->color_info.name.c_str()) == 0 &&
				(current_shader->have_color3 || current_shader->have_color4))
			continue;
		if (strcmp(aname, current_shader->texcoord_info.name.c_str()) == 0 && current_shader->have_texcoord)
			continue;
		eprintf("GLManager::shader_ready : shader %s failed to set attribute %s\n",
			current_shader->name.c_str(),
			current_shader->attribute_info[i].name.c_str());
		ready = false;
	}

	return ready;
}


// Render
void GLManager::draw_points(size_t start, size_t npoints)
{
	// Check shader, but don't refuse to render
	shader_ready();
	glDrawArrays(GL_POINTS, start, npoints);
}

void GLManager::draw_tris(const unsigned *p, std::size_t start, std::size_t ntris)
{
	shader_ready();
	glDrawElements(GL_TRIANGLES, 3 * ntris, GL_UNSIGNED_INT, p + start);
}

void GLManager::draw_tris(std::size_t start, std::size_t ntris)
{
	shader_ready();
	glDrawElements(GL_TRIANGLES, 3 * ntris, GL_UNSIGNED_INT,
		GL_VBO_OFFSET(start * sizeof(unsigned)));
}

void GLManager::draw_tstrips(const unsigned *p, std::size_t start, std::size_t nstripinds)
{
	shader_ready();
	if (enable_restart()) {
		glDrawElements(GL_TRIANGLE_STRIP, nstripinds, GL_UNSIGNED_INT,
			p + start);
	} else {
		p += start;
		const unsigned *end = p + nstripinds;
		while (p < end) {
			const unsigned *stripstart = p;
			while (*p != 0xffffffffu)
				p++;
			glDrawElements(GL_TRIANGLE_STRIP, p - stripstart,
				GL_UNSIGNED_INT, stripstart);
			p++;
		}
	}
}

void GLManager::draw_tstrips(std::size_t start, std::size_t nstripinds)
{
	shader_ready();
	if (enable_restart()) {
		glDrawElements(GL_TRIANGLE_STRIP, nstripinds, GL_UNSIGNED_INT,
			GL_VBO_OFFSET(start * sizeof(unsigned)));
	} else {
		unsigned buf = 0;
		const unsigned *p = 0;
		glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, (GLint *) &buf);
		for (size_t i = 0; i < ibuffers.size(); i++) {
			if (ibuffers[i]->ind == buf) {
				p = (const unsigned *) ibuffers[i]->p;
				break;
			}
		}
		if (!p) {
			eprintf("GLManager::draw_tstrips : drawing tstrips from IBO without restart not supported\n");
			return;
		}
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		draw_tstrips(p, start, nstripinds);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf);
	}
}


// Clear everything we know about
void GLManager::clear()
{
	clear_defines();
	clear_shaders();
	clear_buffers();
	clear_ibuffers();
	clear_textures();
}


// Name lookups
unsigned GLManager::shader(const char *name, bool print_if_not_found)
{
	for (size_t i = 1; i < shaders.size(); i++) {
		if (strcmp(name, shaders[i]->name.c_str()) == 0)
			return shaders[i]->program;
	}
	if (print_if_not_found)
		eprintf("GLManager::shader : lookup of %s failed\n", name);
	return 0;
}

int GLManager::uniform(const char *name, bool print_if_not_found)
{
	if (!current_shader->program) {
		eprintf("GLManager::uniform : no active shader program\n");
		return -1;
	}
	int ret = glGetUniformLocation(current_shader->program, name);
	if (ret < 0 && print_if_not_found)
		eprintf("GLManager::uniform : lookup of %s failed\n", name);
	return ret;
}

int GLManager::attribute(const char *name, bool print_if_not_found)
{
	if (!current_shader->program) {
		eprintf("GLManager::attribute : no active shader program\n");
		return -1;
	}
	int ret = glGetAttribLocation(current_shader->program, name);
	if (ret < 0 && print_if_not_found)
		eprintf("GLManager::attribute : lookup of %s failed\n", name);
	return ret;
}

unsigned GLManager::buffer(const void *p, bool print_if_not_found)
{
	for (size_t i = 1; i < buffers.size(); i++) {
		if (buffers[i]->p == p)
			return buffers[i]->ind;
	}
	if (print_if_not_found)
		eprintf("GLManager::buffer : lookup of %p failed\n", p);
	return 0;
}

unsigned GLManager::ibuffer(const void *p, bool print_if_not_found)
{
	for (size_t i = 1; i < ibuffers.size(); i++) {
		if (ibuffers[i]->p == p)
			return ibuffers[i]->ind;
	}
	if (print_if_not_found)
		eprintf("GLManager::ibuffer : lookup of %p failed\n", p);
	return 0;
}

unsigned GLManager::texture(const void *p, bool print_if_not_found)
{
	for (size_t i = 1; i < textures.size(); i++) {
		if (textures[i]->p == p)
			return textures[i]->ind;
	}
	if (print_if_not_found)
		eprintf("GLManager::texture : lookup of %p failed\n", p);
	return 0;
}


// OpenGL version
float GLManager::glversion()
{
	const char *s = (const char *) glGetString(GL_VERSION);
	const char *c = s;
	while (*c && !isdigit(*c))
		c++;
	float version = 0.0f;
	sscanf(c, "%f", &version);

	return version;
}


// Is this OpenGL ES?
bool GLManager::is_gles()
{
	const char *s = (const char *) glGetString(GL_VERSION);
	return (strncmp(s, "OpenGL ES", 9) == 0);
}


// Is an extension supported?
bool GLManager::extension_supported(const char *name)
{
	const char *ext = (const char *) glGetString(GL_EXTENSIONS);
	const char *ext_end = ext + strlen(ext);
	size_t namelen = strlen(name);

	for (;;) {
		ext = strstr(ext, name);
		if (!ext || ext + namelen > ext_end)
			break;
		if (*(ext + namelen) == ' ' || *(ext + namelen) == '\0')
			return true;
		ext++;
	}
	return false;
}


// Do we have primitive restart?
bool GLManager::have_primitive_restart()
{
	check_restart();
	return (have_restart_fixed || have_restart || have_restart_nv);
}


// Is rendering of tstrips slow?
bool GLManager::slow_tstrips()
{
	// Intel GPUs without primitive restart
	return strstr((const char *) glGetString(GL_RENDERER), "Intel") &&
	       !have_primitive_restart();
}


// Check whether we have primitive restart in some form
void GLManager::check_restart()
{
	if (checked_restart)
		return;
	checked_restart = true;

	float version = glversion();
	bool gles = is_gles();
	have_restart_fixed = ((gles && version >= 3.0f) || version >= 4.3f);
	have_restart = (!gles && version >= 3.1f);
	have_restart_nv = (extension_supported("GL_NV_primitive_restart") &&
		glGetProcs() && glPrimitiveRestartIndexNV != NULL);
}


// Turn on primitive restart
bool GLManager::enable_restart()
{
	check_restart();
	if (have_restart_fixed) {
		glEnable(GL_PRIMITIVE_RESTART_FIXED_INDEX);
		return true;
	} else if (have_restart) {
		glPrimitiveRestartIndex(0xffffffffu);
		glEnable(GL_PRIMITIVE_RESTART);
		return true;
	} else if (have_restart_nv) {
		glPrimitiveRestartIndexNV(0xffffffffu);
		glEnableClientState(GL_PRIMITIVE_RESTART_NV);
		return true;
	} else {
		return false;
	}
}


// Check whether we have glClampColor or glClampColorARB, and if so
// set clamping off
void GLManager::check_clamp()
{
	if (checked_clamp)
		return;
	checked_clamp = true;

	float version = glversion();
	if (version >= 3.0f) {
		glGetProcs();
		glClampColor(GL_CLAMP_VERTEX_COLOR, GL_FALSE);
		glClampColor(GL_CLAMP_FRAGMENT_COLOR, GL_FALSE);
		glClampColor(GL_CLAMP_READ_COLOR, GL_FALSE);
	} else if (extension_supported("GL_ARB_color_buffer_float") &&
	           glGetProcs() && glClampColorARB != NULL) {
		glClampColorARB(GL_CLAMP_VERTEX_COLOR, GL_FALSE);
		glClampColorARB(GL_CLAMP_FRAGMENT_COLOR, GL_FALSE);
		glClampColorARB(GL_CLAMP_READ_COLOR, GL_FALSE);
	}
}


// Check whether we have Vertex Array Objects, and if so set one up
void GLManager::check_vao()
{
	if (checked_vao)
		return;
	checked_vao = true;

	float version = glversion();
	if (version >= 3.0f) {  // Both GL and GLES
		glGetProcs();
		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);
	}
}


// Check whether a shader compiled OK
bool GLManager::shader_ok(unsigned shader_num, const char *what, const char *name)
{
	GLint status;
	glGetShaderiv(shader_num, GL_COMPILE_STATUS, &status);
	if (status == GL_TRUE)
		return true;

	GLint log_length;
	glGetShaderiv(shader_num, GL_INFO_LOG_LENGTH, &log_length);
	vector<char> log(log_length);
	glGetShaderInfoLog(shader_num, log_length, NULL, &log[0]);

	eprintf("GLManager::make_shader : Error %s %s\n%s\n",
		what, name, &log[0]);
	return false;
}


// Check whether a shader program linked OK
bool GLManager::program_ok(unsigned prog_num, const char *what, const char *name)
{
	GLint status;
	glGetProgramiv(prog_num, GL_LINK_STATUS, &status);

	GLint log_length;
	glGetProgramiv(prog_num, GL_INFO_LOG_LENGTH, &log_length);
	if (log_length > 5) {
		vector<char> log(log_length);
		glGetProgramInfoLog(prog_num, log_length, NULL, &log[0]);

		eprintf("GLManager::make_shader : %s %s %s\n%s\n",
			status ? "Warnings" : "Error", what, name, &log[0]);
	}
	return status;
}


// Set up the hook for eprintf
void GLManager::set_error_hook(void (*hook)(const char *))
{
	eprintf_hook = hook;
}


// Error printout
void GLManager::eprintf(const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	char buf[1024];
	vsnprintf(buf, sizeof(buf), format, ap);
	va_end(ap);

	if (eprintf_hook) {
		eprintf_hook(buf);
	} else {
		fputs(buf, stderr);
		fflush(stderr);
	}
}

} // namespace trimesh
