#ifndef SHADERS_INC
#define SHADERS_INC
/*
Szymon Rusinkiewicz
Princeton University

shaders.inc
GLSL shaders for mesh_view
*/


#ifndef GL_VERTEX_PROGRAM_POINT_SIZE
#define GL_VERTEX_PROGRAM_POINT_SIZE 0x8642
#endif


// Unlit - just pass through gl_Color
const char unlit_vert[] = shader_text(
	varying vec4 color;
	void main()
	{
		gl_Position = ftransform();
		color = gl_Color;
	}
);

const char unlit_point_vert[] = shader_text(
	uniform float pointscale;
	varying vec4 color;
	void main()
	{
		gl_Position = ftransform();
		vec4 worldPos = gl_ModelViewMatrix * gl_Vertex;
		float z = -worldPos.z / worldPos.w;
		gl_PointSize = max(1.0, pointscale / z);
		color = gl_Color;
	}
);

const char unlit_frag[] = shader_text(
	varying vec4 color;
	void main()
	{
		gl_FragColor = color;
	}
);


// Lit with per-fragment Phong based on interpolated normals
const char phong_vert[] = shader_text(
	varying vec3 norm;
	varying vec4 color;
	void main()
	{
		gl_Position = ftransform();
		norm = (gl_NormalMatrix * gl_Normal).xyz;
		color = gl_Color;
	}
);

const char phong_point_vert[] = shader_text(
	uniform float pointscale;
	varying vec3 norm;
	varying vec4 color;
	void main()
	{
		gl_Position = ftransform();
		vec4 worldPos = gl_ModelViewMatrix * gl_Vertex;
		float z = -worldPos.z / worldPos.w;
		gl_PointSize = max(1.0, pointscale / z);
		norm = (gl_NormalMatrix * gl_Normal).xyz;
		color = gl_Color;
	}
);

const char phong_frag[] = shader_text(
	varying vec3 norm;
	varying vec4 color;
	void main()
	{
		vec4 ambient = color * gl_LightSource[0].ambient;
		vec3 nnorm = normalize(norm);
		float ndotl = abs(dot(nnorm, gl_LightSource[0].position.xyz));
		vec4 diffuse = color * gl_LightSource[0].diffuse * ndotl;
		float ndoth = abs(dot(nnorm, gl_LightSource[0].halfVector.xyz));
		vec4 specular = gl_FrontMaterial.specular *
			pow(ndoth, gl_FrontMaterial.shininess);
		gl_FragColor = ambient + diffuse + specular;
	}
);


// Lit with per-fragment Phong based on (flat) face normals
const char flat_vert[] = shader_text(
	varying vec3 pos;
	varying vec4 color;
	void main()
	{
		gl_Position = ftransform();
		pos = (gl_ModelViewMatrix * gl_Vertex).xyz;
		color = gl_Color;
	}
);

const char flat_frag[] = shader_text(
	varying vec3 pos;
	varying vec4 color;
	void main()
	{
		vec4 ambient = color * gl_LightSource[0].ambient;
		vec3 dposdx = vec3(dFdx(pos.x), dFdx(pos.y), dFdx(pos.z));
		vec3 dposdy = vec3(dFdy(pos.x), dFdy(pos.y), dFdy(pos.z));
		vec3 nnorm = normalize(cross(dposdx, dposdy));
		float ndotl = abs(dot(nnorm, gl_LightSource[0].position.xyz));
		vec4 diffuse = color * gl_LightSource[0].diffuse * ndotl;
		float ndoth = abs(dot(nnorm, gl_LightSource[0].halfVector.xyz));
		vec4 specular = gl_FrontMaterial.specular *
			pow(ndoth, gl_FrontMaterial.shininess);
		gl_FragColor = ambient + diffuse + specular;
	}
);

#endif
