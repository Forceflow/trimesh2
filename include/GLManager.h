#ifndef GLMANAGER_H
#define GLMANAGER_H
/*
Szymon Rusinkiewicz
Princeton University

GLManager.h
Manager for modern OpenGL buffers and shaders.
*/

#include <cstddef>
#include <vector>

namespace trimesh {

using ::std::size_t;

class GLManager {
private:
	bool glGetProcs();

public:
	GLManager();
	~GLManager();

	// #define used by the next make_shader
	void define(const char *defname, const char *value);

	// Clear all defines.  Note that make_shader() does *not* do this.
	void clear_defines();

	// Create a shader program using the given vertex and fragment shaders
	unsigned make_shader(const char *name,
		const char *vert_text, const char *frag_text);
	unsigned make_shader_from_files(const char *name,
		const char *vert_filename, const char *frag_filename);

	// Include shader code directly in source file, via stringification.
	// Note that line numbering for errors will not work.  Usage:
	// const char my_vert_shader[] = shader_text(
	// 	void main()
	//	{
	//		gl_Position = ftransform();
	//	}
	// );
#define shader_text(shader) #shader

	// Make shader active, for rendering and uniform/attribute definitions
	bool use_shader(unsigned ind = 0);
	bool use_shader(const char *name)
		{ return use_shader(shader(name, true)); }

	// Delete a shader, or all shaders
	void clear_shader(unsigned ind);
	void clear_shader(const char *name)
		{ clear_shader(shader(name, true)); }
	void clear_shaders();

	// Set a uniform for the current shader, by location or name
private:
	bool uniform_boilerplate(int ind, unsigned type);
public:
	void uniform1i(int ind, int x);
	void uniform2i(int ind, int x, int y);
	void uniform3i(int ind, int x, int y, int z);
	void uniform4i(int ind, int x, int y, int z, int w);
	void uniform1iv(int ind, const int *p)
		{ uniform1i(ind, p[0]); }
	void uniform2iv(int ind, const int *p)
		{ uniform2i(ind, p[0], p[1]); }
	void uniform3iv(int ind, const int *p)
		{ uniform3i(ind, p[0], p[1], p[2]); }
	void uniform4iv(int ind, const int *p)
		{ uniform4i(ind, p[0], p[1], p[2], p[3]); }
	void uniform1f(int ind, float x);
	void uniform2f(int ind, float x, float y);
	void uniform3f(int ind, float x, float y, float z);
	void uniform4f(int ind, float x, float y, float z, float w);
	void uniform1fv(int ind, const float *p)
		{ uniform1f(ind, p[0]); }
	void uniform2fv(int ind, const float *p)
		{ uniform2f(ind, p[0], p[1]); }
	void uniform3fv(int ind, const float *p)
		{ uniform3f(ind, p[0], p[1], p[2]); }
	void uniform4fv(int ind, const float *p)
		{ uniform4f(ind, p[0], p[1], p[2], p[3]); }
	void uniformMatrix2fv(int ind, const float *p);
	void uniformMatrix3fv(int ind, const float *p);
	void uniformMatrix4fv(int ind, const float *p);
	void uniformTexture(int ind, unsigned texind);
	void uniformTexture(int ind, const void *p)
		{ uniformTexture(ind, texture(p, true)); }
	template <class T> void uniformTexture(int ind, const ::std::vector<T> &v)
		{ uniformTexture(ind, texture(v, true)); }

	void uniform1i(const char *name, int x)
		{ uniform1i(uniform(name, true), x); }
	void uniform2i(const char *name, int x, int y)
		{ uniform2i(uniform(name, true), x, y); }
	void uniform3i(const char *name, int x, int y, int z)
		{ uniform3i(uniform(name, true), x, y, z); }
	void uniform4i(const char *name, int x, int y, int z, int w)
		{ uniform4i(uniform(name, true), x, y, z, w); }
	void uniform1iv(const char *name, const int *p)
		{ uniform1iv(uniform(name, true), p); }
	void uniform2iv(const char *name, const int *p)
		{ uniform2iv(uniform(name, true), p); }
	void uniform3iv(const char *name, const int *p)
		{ uniform3iv(uniform(name, true), p); }
	void uniform4iv(const char *name, const int *p)
		{ uniform4iv(uniform(name, true), p); }
	void uniform1f(const char *name, float x)
		{ uniform1f(uniform(name, true), x); }
	void uniform2f(const char *name, float x, float y)
		{ uniform2f(uniform(name, true), x, y); }
	void uniform3f(const char *name, float x, float y, float z)
		{ uniform3f(uniform(name, true), x, y, z); }
	void uniform4f(const char *name, float x, float y, float z, float w)
		{ uniform4f(uniform(name, true), x, y, z, w); }
	void uniform1fv(const char *name, const float *p)
		{ uniform1fv(uniform(name, true), p); }
	void uniform2fv(const char *name, const float *p)
		{ uniform2fv(uniform(name, true), p); }
	void uniform3fv(const char *name, const float *p)
		{ uniform3fv(uniform(name, true), p); }
	void uniform4fv(const char *name, const float *p)
		{ uniform4fv(uniform(name, true), p); }
	void uniformMatrix2fv(const char *name, const float *p)
		{ uniformMatrix2fv(uniform(name, true), p); }
	void uniformMatrix3fv(const char *name, const float *p)
		{ uniformMatrix3fv(uniform(name, true), p); }
	void uniformMatrix4fv(const char *name, const float *p)
		{ uniformMatrix4fv(uniform(name, true), p); }
	void uniformTexture(const char *name, unsigned texind)
		{ uniformTexture(uniform(name, true), texind); }
	void uniformTexture(const char *name, const void *p)
		{ uniformTexture(uniform(name, true), texture(p, true)); }
	template <class T> void uniformTexture(const char *name, const ::std::vector<T> &v)
		{ uniformTexture(uniform(name, true), texture(v, true)); }

	// "Clear" uniforms so that they have to be set again
	void clear_uniform(int ind);
	void clear_uniform(const char *name)
		{ clear_uniform(uniform(name, true)); }
	void clear_uniforms();

	// Set an attribute to a (client-side) array or vector
private:
	void attributearray(int ind, int floats_per_attribute, const float *p,
		size_t stride = 0, size_t offset = 0);
public:
	void attributearray1fv(int ind, const float *p,
		size_t stride = 0, size_t offset = 0)
		{ attributearray(ind, 1, p, stride, offset); }
	void attributearray2fv(int ind, const float *p,
		size_t stride = 0, size_t offset = 0)
		{ attributearray(ind, 2, p, stride, offset); }
	void attributearray3fv(int ind, const float *p,
		size_t stride = 0, size_t offset = 0)
		{ attributearray(ind, 3, p, stride, offset); }
	void attributearray4fv(int ind, const float *p,
		size_t stride = 0, size_t offset = 0)
		{ attributearray(ind, 4, p, stride, offset); }

	void attributearray1fv(const char *name, const float *p,
		size_t stride = 0, size_t offset = 0)
		{ attributearray1fv(attribute(name, true), p, stride, offset); }
	void attributearray2fv(const char *name, const float *p,
		size_t stride = 0, size_t offset = 0)
		{ attributearray2fv(attribute(name, true), p, stride, offset); }
	void attributearray3fv(const char *name, const float *p,
		size_t stride = 0, size_t offset = 0)
		{ attributearray3fv(attribute(name, true), p, stride, offset); }
	void attributearray4fv(const char *name, const float *p,
		size_t stride = 0, size_t offset = 0)
		{ attributearray4fv(attribute(name, true), p, stride, offset); }

	template <class T> void attributearray1fv(int ind, const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ attributearray1fv(ind, (const float *) &v[0], stride, offset); }
	template <class T> void attributearray2fv(int ind, const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ attributearray2fv(ind, (const float *) &v[0], stride, offset); }
	template <class T> void attributearray3fv(int ind, const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ attributearray3fv(ind, (const float *) &v[0], stride, offset); }
	template <class T> void attributearray4fv(int ind, const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ attributearray4fv(ind, (const float *) &v[0], stride, offset); }

	template <class T> void attributearray1fv(const char *name, const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ attributearray1fv(attribute(name, true), (const float *) &v[0], stride, offset); }
	template <class T> void attributearray2fv(const char *name, const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ attributearray2fv(attribute(name, true), (const float *) &v[0], stride, offset); }
	template <class T> void attributearray3fv(const char *name, const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ attributearray3fv(attribute(name, true), (const float *) &v[0], stride, offset); }
	template <class T> void attributearray4fv(const char *name, const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ attributearray4fv(attribute(name, true), (const float *) &v[0], stride, offset); }

	// Set an attribute to a buffer
private:
	void attributearray(int ind, int floats_per_attribute, unsigned buf,
		size_t stride = 0, size_t offset = 0);
public:
	void attributearray1fv(int ind, unsigned buf,
		size_t stride = 0, size_t offset = 0)
		{ attributearray(ind, 1, buf, stride, offset); }
	void attributearray2fv(int ind, unsigned buf,
		size_t stride = 0, size_t offset = 0)
		{ attributearray(ind, 2, buf, stride, offset); }
	void attributearray3fv(int ind, unsigned buf,
		size_t stride = 0, size_t offset = 0)
		{ attributearray(ind, 3, buf, stride, offset); }
	void attributearray4fv(int ind, unsigned buf,
		size_t stride = 0, size_t offset = 0)
		{ attributearray(ind, 4, buf, stride, offset); }

	void attributearray1fv(const char *name, unsigned buf,
		size_t stride = 0, size_t offset = 0)
		{ attributearray1fv(attribute(name, true), buf, stride, offset); }
	void attributearray2fv(const char *name, unsigned buf,
		size_t stride = 0, size_t offset = 0)
		{ attributearray2fv(attribute(name, true), buf, stride, offset); }
	void attributearray3fv(const char *name, unsigned buf,
		size_t stride = 0, size_t offset = 0)
		{ attributearray3fv(attribute(name, true), buf, stride, offset); }
	void attributearray4fv(const char *name, unsigned buf,
		size_t stride = 0, size_t offset = 0)
		{ attributearray4fv(attribute(name, true), buf, stride, offset); }

	// Set an attribute to a constant
private:
	bool attribute_boilerplate(int ind, int floats_per_attribute);
public:
	void attribute1f(int ind, float x);
	void attribute2f(int ind, float x, float y);
	void attribute3f(int ind, float x, float y, float z);
	void attribute4f(int ind, float x, float y, float z, float w);
	void attribute1fv(int ind, const float *p)
		{ uniform1f(ind, p[0]); }
	void attribute2fv(int ind, const float *p)
		{ uniform2f(ind, p[0], p[1]); }
	void attribute3fv(int ind, const float *p)
		{ uniform3f(ind, p[0], p[1], p[2]); }
	void attribute4fv(int ind, const float *p)
		{ uniform4f(ind, p[0], p[1], p[2], p[3]); }
	void attribute1f(const char *name, float x)
		{ attribute1f(attribute(name, true), x); }
	void attribute2f(const char *name, float x, float y)
		{ attribute2f(attribute(name, true), x, y); }
	void attribute3f(const char *name, float x, float y, float z)
		{ attribute3f(attribute(name, true), x, y, z); }
	void attribute4f(const char *name, float x, float y, float z, float w)
		{ attribute4f(attribute(name, true), x, y, z, w); }
	void attribute1fv(const char *name, const float *p)
		{ attribute1fv(attribute(name, true), p); }
	void attribute2fv(const char *name, const float *p)
		{ attribute2fv(attribute(name, true), p); }
	void attribute3fv(const char *name, const float *p)
		{ attribute3fv(attribute(name, true), p); }
	void attribute4fv(const char *name, const float *p)
		{ attribute4fv(attribute(name, true), p); }

	// Clear attribute bindings
	void clear_attribute(int ind);
	void clear_attribute(const char *name)
		{ clear_attribute(attribute(name, true)); }
	void clear_attributes();

	// Old-style attributes - in pointer, vector, and buffer variants
	void vertexarray3fv(const float *p,
		size_t stride = 0, size_t offset = 0);
	void normalarray3fv(const float *p,
		size_t stride = 0, size_t offset = 0);
	void colorarray3fv(const float *p,
		size_t stride = 0, size_t offset = 0);
	void colorarray4fv(const float *p,
		size_t stride = 0, size_t offset = 0);
	void texcoordarray2fv(const float *p,
		size_t stride = 0, size_t offset = 0);

	template <class T> void vertexarray3fv(const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ vertexarray3fv((const float *) &v[0], stride, offset); }
	template <class T> void normalarray3fv(const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ normalarray3fv((const float *) &v[0], stride, offset); }
	template <class T> void colorarray3fv(const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ colorarray3fv((const float *) &v[0], stride, offset); }
	template <class T> void colorarray4fv(const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ colorarray4fv((const float *) &v[0], stride, offset); }
	template <class T> void texcoordarray2fv(const ::std::vector<T> &v,
		size_t stride = 0, size_t offset = 0)
		{ texcoordarray2fv((const float *) &v[0], stride, offset); }

	void vertexarray3fv(unsigned buf,
		size_t stride = 0, size_t offset = 0);
	void normalarray3fv(unsigned buf,
		size_t stride = 0, size_t offset = 0);
	void colorarray3fv(unsigned buf,
		size_t stride = 0, size_t offset = 0);
	void colorarray4fv(unsigned buf,
		size_t stride = 0, size_t offset = 0);
	void texcoordarray2fv(unsigned buf,
		size_t stride = 0, size_t offset = 0);

	void clear_vertexarray();
	void clear_normalarray();
	void clear_colorarray();
	void clear_texcoordarray();

	// Old-style constant attributes
	void normal3f(float nx, float ny, float nz);
	void normal3fv(const float *p)
		{ normal3f(p[0], p[1], p[2]); }
	void color3f(float r, float g, float b);
	void color3fv(const float *p)
		{ color3f(p[0], p[1], p[2]); }
	void color4f(float r, float g, float b, float a);
	void color4fv(const float *p)
		{ color4f(p[0], p[1], p[2], p[3]); }
	void texcoord2f(float u, float v);
	void texcoord2fv(const float *p)
		{ texcoord2f(p[0], p[1]); }

	// Array buffers
	unsigned make_buffer(const float *p, size_t nfloat);
	template <class T> unsigned make_buffer(const ::std::vector<T> &v)
		{ return make_buffer((const float *) &v[0], v.size() * sizeof(v[0]) / sizeof(float)); }
	void update_buffer(unsigned buf, const float *p, size_t nfloat);
	template <class T> void update_buffer(unsigned buf, const ::std::vector<T> &v)
		{ update_buffer(buf, (const float *) &v[0], v.size() * sizeof(v[0]) / sizeof(float)); }
	bool use_buffer(unsigned buf = 0);
	bool use_buffer(const void *p)
		{ return use_buffer(buffer(p)); }
	template <class T> bool use_buffer(const ::std::vector<T> &v)
		{ return use_buffer(buffer(v)); }
	void clear_buffer(unsigned buf);
	void clear_buffer(const void *p)
		{ clear_buffer(buffer(p)); }
	template <class T> void clear_buffer(const ::std::vector<T> &v)
		{ clear_buffer(buffer(v)); }
	void clear_buffers();

	// Index buffers
	unsigned make_ibuffer(const unsigned *p, size_t ninds);
	unsigned make_ibuffer(const int *p, size_t ninds)
		{ return make_ibuffer((const unsigned *) p, ninds); }
	template <class T> unsigned make_ibuffer(const ::std::vector<T> &v)
		{ return make_ibuffer((const unsigned *) &v[0], v.size() * sizeof(v[0]) / sizeof(unsigned)); }
	void update_ibuffer(unsigned buf, const unsigned *p, size_t ninds);
	void update_ibuffer(unsigned buf, const int *p, size_t ninds)
		{ update_ibuffer(buf, (const unsigned *) p, ninds); }
	template <class T> void update_ibuffer(unsigned buf, const ::std::vector<T> &v)
		{ update_ibuffer(buf, (const unsigned *) &v[0], v.size() * sizeof(v[0]) / sizeof(unsigned)); }
	bool use_ibuffer(unsigned buf = 0);
	bool use_ibuffer(const void *p)
		{ return use_ibuffer(ibuffer(p)); }
	template <class T> bool use_ibuffer(const ::std::vector<T> &v)
		{ return use_ibuffer(ibuffer(v)); }
	void clear_ibuffer(unsigned buf);
	void clear_ibuffer(const void *p)
		{ clear_ibuffer(buffer(p)); }
	template <class T> void clear_ibuffer(const ::std::vector<T> &v)
		{ clear_ibuffer(buffer(v)); }
	void clear_ibuffers();

private:
	// Length of a pixel in bytes
	int pixel_in_bytes(unsigned format, unsigned type);

public:
	// Should future texture loads flip the image in Y?
	void flip_textures(bool flip);

	// Textures
	unsigned make_texture(int w, int h, const void *p,
		unsigned format = 0x1907u /* GL_RGB */,
		unsigned type = 0x1401u /* GL_UNSIGNED_BYTE */,
		unsigned internalformat = 0 /* defaults to format */);
	template <class T> unsigned make_texture(int w, int h,
		const ::std::vector<T> &v, unsigned format = 0x1907u,
		unsigned type = 0x1401u, unsigned internalformat = 0)
		{ return make_texture(w, h, (const void *) &v[0], format, type, internalformat); }
	void update_texture(unsigned texind, int w, int h, const void *p,
		unsigned format = 0x1907u, unsigned type = 0x1401u,
		unsigned internalformat = 0);
	template <class T> void update_texture(unsigned texind, int w, int h,
		const ::std::vector<T> &v, unsigned format = 0x1907u,
		unsigned type = 0x1401u, unsigned internalformat = 0)
		{ update_texture(texind, w, h, (const void *) &v[0], format, type, internalformat); }

	bool use_texture(unsigned texind = 0);
	bool use_texture(const void *p)
		{ return use_texture(texture(p, true)); }
	template <class T> bool use_texture(const ::std::vector<T> &v)
		{ return use_texture(texture(v, true)); }

	void clear_texture(unsigned ind);
	void clear_texture(const void *p)
		{ clear_texture(texture(p, true)); }
	template <class T> void clear_texture(const ::std::vector<T> &v)
		{ clear_texture(texture(v, true)); }
	void clear_textures();

	// Render to texture
	void rendertarget(unsigned texind = 0, unsigned depthbufferind = 0);

	// Is shader program ready to go, and have all required
	// attributes and uniforms been specified?
	bool shader_ready();

	// Render
	void draw_points(size_t start, size_t npoints);
	void draw_points(size_t npoints)
		{ draw_points(0, npoints); }

	void draw_tris(const unsigned *p, size_t start, size_t ntris);
	void draw_tris(const unsigned *p, size_t ntris)
		{ draw_tris(p, 0, ntris); }
	void draw_tris(const int *p, size_t start, size_t ntris)
		{ draw_tris((const unsigned *) p, start, ntris); }
	void draw_tris(const int *p, size_t ntris)
		{ draw_tris((const unsigned *) p, 0, ntris); }
	template <class T> void draw_tris(const ::std::vector<T> &v)
		{ draw_tris((const unsigned *) &v[0], 0, v.size()); }
	template <class T> void draw_tris(const ::std::vector<T> &v, size_t start, size_t ntris)
		{ draw_tris((const unsigned *) &v[0], start, ntris); }
	void draw_tris(size_t start, size_t ntris);
	void draw_tris(int start, size_t ntris)
		{ draw_tris((size_t) start, ntris); }
	void draw_tris(size_t ntris)
		{ draw_tris((size_t) 0, ntris); }

	void draw_tstrips(const unsigned *p, size_t start, size_t nstripinds);
	void draw_tstrips(const unsigned *p, size_t nstripinds)
		{ draw_tstrips(p, 0, nstripinds); }
	void draw_tstrips(const int *p, size_t start, size_t nstripinds)
		{ draw_tstrips((const unsigned *) p, start, nstripinds); }
	void draw_tstrips(const int *p, size_t nstripinds)
		{ draw_tstrips((const unsigned *) p, 0, nstripinds); }
	template <class T> void draw_tstrips(const ::std::vector<T> &v)
		{ draw_tstrips((const unsigned *) &v[0], 0, v.size()); }
	template <class T> void draw_tstrips(const ::std::vector<T> &v, size_t start, size_t nstripinds)
		{ draw_tstrips((const unsigned *) &v[0], start, nstripinds); }
	void draw_tstrips(size_t start, size_t nstripinds);
	void draw_tstrips(int start, size_t nstripinds)
		{ draw_tstrips((size_t) start, nstripinds); }
	void draw_tstrips(size_t nstripinds)
		{ draw_tstrips((size_t) 0, nstripinds); }

	// Clear everything we know about
	void clear();

	// Name lookups
	unsigned shader(const char *name, bool print_if_not_found = false);
	int uniform(const char *name, bool print_if_not_found = false);
	int attribute(const char *name, bool print_if_not_found = false);
	unsigned buffer(const void *p, bool print_if_not_found = false);
	template <class T> unsigned buffer(const ::std::vector<T> &v, bool print_if_not_found = false)
		{ return buffer((const void *) &v[0], print_if_not_found); }
	unsigned ibuffer(const void *p, bool print_if_not_found = false);
	template <class T> unsigned ibuffer(const ::std::vector<T> &v, bool print_if_not_found = false)
		{ return ibuffer((const void *) &v[0], print_if_not_found); }
	unsigned texture(const void *p, bool print_if_not_found = false);
	template <class T> unsigned texture(const ::std::vector<T> &v, bool print_if_not_found = false)
		{ return texture((const void *) &v[0], print_if_not_found); }

	// OpenGL version
	float glversion();

	// Is this OpenGL ES?
	bool is_gles();

	// Is an extension supported?
	bool extension_supported(const char *name);

	// Check for primitive restart
	bool have_primitive_restart();

	// Is rendering of tstrips slow?
	bool slow_tstrips();

	// Set hook for error printouts
	void set_error_hook(void (*hook)(const char *));

private:
	class ShaderDefine;
	struct UniformInfo;
	struct AttributeInfo;
	class ShaderInfo;
	class BufferInfo;
	::std::vector<ShaderDefine *> defines;
	::std::vector<ShaderInfo *> shaders;
	::std::vector<BufferInfo *> buffers, ibuffers, textures;
	ShaderInfo *current_shader;
	BufferInfo *current_buffer, *current_ibuffer, *current_texture;
	bool flipping_textures;
	bool checked_restart, have_restart_fixed, have_restart, have_restart_nv;
	void check_restart();
	bool enable_restart();
	bool checked_clamp;
	void check_clamp();
	bool checked_vao;
	unsigned vao;
	void check_vao();
	bool shader_ok(unsigned shader_num, const char *what, const char *name);
	bool program_ok(unsigned prog_num, const char *what, const char *name);
	void eprintf(const char *format, ...);
	void (*eprintf_hook)(const char *);
};

} // namespace trimesh

#endif
