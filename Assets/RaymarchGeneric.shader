﻿Shader "Hidden/RaymarchGeneric"
{
		SubShader
	{
		// No culling or depth
		Cull Off ZWrite Off ZTest Always

		Pass
		{
			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag
			// You may need to use an even later shader model target, depending on how many instructions you have
			// or if you need variable-length for loops.
			#pragma target 3.0

			#include "UnityCG.cginc"
			#include "DistanceFunc.cginc"
			
			uniform sampler2D _CameraDepthTexture;
			uniform float4x4 _CameraClipToWorld;
			// These three are set by our script (see RaymarchGeneric.cs)
			uniform sampler2D _MainTex;
			uniform float4x4 _FrustumCornersWS;
			uniform float4 _CameraWS;
			uniform float3 _LightDir;

			struct appdata
			{
				// Remember, the z value here contains the index of _FrustumCornersWS to use
				float4 vertex : POSITION;
				float2 uv : TEXCOORD0;
			};

			struct v2f
			{
				float4 pos : SV_POSITION;
				float2 uv : TEXCOORD0;
				float4 ray : TEXCOORD1;
			};

			v2f vert (appdata v)
			{
				v2f o;
				
				half index = v.vertex.z;
				v.vertex.z = 0.1;
				
				o.pos = mul(UNITY_MATRIX_MVP, v.vertex);
				o.uv = v.uv.xy;
				
				#if UNITY_UV_STARTS_AT_TOP
				o.uv.y = 1 - o.uv.y;
				#endif

				o.ray = _FrustumCornersWS[(int)index];

				return o;
			}

			float map(float3 p) {
				return sdTorus(p - float3(_Time.y,0,0), float2(1, 0.2));
			}

			float3 calcNormal(in float3 pos)
			{
				const float2 eps = float2(0.001, 0.0);
				// The idea here is to find the "gradient" of the distance field at pos
				// Remember, the distance field is not boolean - even if you are inside an object
				// the number is negative, so this calculation still works.
				// Essentially you are approximating the derivative of the distance field at this point.
				float3 nor = float3(
					map(pos + eps.xyy).x - map(pos - eps.xyy).x,
					map(pos + eps.yxy).x - map(pos - eps.yxy).x,
					map(pos + eps.yyx).x - map(pos - eps.yyx).x);
				return normalize(nor);
			}

			// Raymarch along given ray
			// ro: ray origin
			// rd: ray direction
			// s: unity depth buffer
			fixed4 raymarch(float3 ro, float3 rd, float s) {
				fixed4 ret = fixed4(0,0,0,0);

				const int maxstep = 64;
				float t = 0; // current distance traveled along ray
				for (int i = 0; i < maxstep; ++i) {
					if (t > s) {
						ret = fixed4(0, 0, 0, 0);
						break;
					}
					float3 p = ro + rd * t;
					float d = map(p);

					if (d < 0.001) {
						float3 n = calcNormal(p);
						ret = fixed4(dot(-_LightDir.xyz, n).rrr, 1);
						break;
					}

					t += min(d, s - t);
				}

				return ret;
			}

			fixed4 frag (v2f i) : SV_Target
			{
				// ray direction
				float3 rd = normalize(i.ray.xyz);
				// ray origin (camera position)
				float3 ro = _CameraWS;

				float2 duv = i.uv;
				#if UNITY_UV_STARTS_AT_TOP
				duv.y = 1 - duv.y;
				#endif

				// Convert from depth buffer to true distance from camera
				float depth = tex2D(_CameraDepthTexture, duv).r;
				float4 projPos = float4(duv.x * 2 - 1, duv.y * 2 - 1, depth * 2 - 1, 1.0f);
				float4 posvs = mul(_CameraClipToWorld, projPos);
				posvs /= posvs.w;
				depth = length(posvs.xyz - ro);

				fixed3 col = tex2D(_MainTex,i.uv);
				fixed4 add = raymarch(ro, rd, depth);

				// Returns final color using alpha blending
				return fixed4(col*(1.0 - add.w) + add.xyz * add.w,1.0);
			}
			ENDCG
		}
	}
}