Shader "GLSL shader with single texture" {
   Properties {
      _MainTex ("Texture Image", 2D) = "white" {} 
   }
   SubShader {
      Pass {	
         GLSLPROGRAM
		 
		float smin(float a, float b, float k)
		{
			float h = clamp(0.5 + 0.5*(b - a) / k, 0.0, 1.0);
			return mix(b, a, h) - k * h*(1.0 - h);
		}

		vec2 smin(vec2 a, vec2 b, float k)
		{
			float h = clamp(0.5 + 0.5*(b.x - a.x) / k, 0.0, 1.0);
			return vec2(mix(b.x, a.x, h) - k * h*(1.0 - h), mix(b.y, a.y, h));
		}
		vec4 smin(vec4 a, vec4 b, float k)
		{
			float h = clamp(0.5 + 0.5*(b.x - a.x) / k, 0.0, 1.0);
			return vec4(mix(b.x, a.x, h) - k * h*(1.0 - h), mix(b.yzw, a.yzw, h));
		}

		float smax(float a, float b, float k)
		{
			float h = clamp(0.5 + 0.5*(b - a) / k, 0.0, 1.0);
			return mix(a, b, h) + k * h*(1.0 - h);
		}

		vec3 smax(vec3 a, vec3 b, float k)
		{
			vec3 h = clamp(0.5 + 0.5*(b - a) / k, 0.0, 1.0);
			return mix(a, b, h) + k * h*(1.0 - h);
		}
                  
         uniform sampler2D _MainTex;	
         uniform vec4 _MainTex_ST; 
            // tiling and offset parameters of property             

         varying vec4 textureCoordinates; 

         #ifdef VERTEX
                  
         void main()
         {
            textureCoordinates = gl_MultiTexCoord0;
            gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
         }
         
         #endif

         #ifdef FRAGMENT
         
         void main()
         {
            gl_FragColor = texture2D(_MainTex, 
               _MainTex_ST.xy * textureCoordinates.xy 
               + _MainTex_ST.zw);	
               // textureCoordinates are multiplied with the tiling 
               // parameters and the offset parameters are added
         }
         
         #endif

         ENDGLSL
      }
   }
   // The definition of a fallback shader should be commented out 
   // during development:
   // Fallback "Unlit/Texture"
}