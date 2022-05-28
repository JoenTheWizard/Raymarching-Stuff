///<music id="PyR8su4x8m4" loop="true">
float D(vec3 p){p+=vec3(0,T,cos(15.*T)+5.*T);p=abs(mod(p,2.)-.5*2.);return(p.x+p.y+p.z-.2*cos(15.*T)-.35)*.4;}void main(){vec3 H=normalize(vec3(L+L-(K.xy=R.xy),K.y));K-=K;for(int i;i++<50;)K+=D(H*K.x);K.rgb=floor(K.rgb/vec3(3)*1.);K=1./K*vec4(cos(T),sin(T),.3,1);}
