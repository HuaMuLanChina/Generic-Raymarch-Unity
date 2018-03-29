Shader "GLSL shader with single texture" {
   Properties {
      _MainTex ("Texture Image", 2D) = "white" {} 
   }
   SubShader {
      Pass {	
         GLSLPROGRAM
                  
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