///<music id="SZN3As_oZwI" loop="true">
float D(vec3 p){
	p.z-=3.;
    float dotp = dot(p,p);
    p.x+=sin(T*40.)*.001;
    p=p/dotp*4.;
    p=sin(p+vec3(sin(1.+T)*2.,-T,-T*2.));
    float d=length(p.yz)-.22;
    d=min(min(min(d,length(p.xz)-.22),length(p.xy)-.22),length(p*p*p)-.22*.3);return d*dotp/6.;}
void main(){vec3 H=normalize(vec3(L+L-(K.xy=R.xy),K.y));K-=K;for(int i;i++<25;)K+=D(H*K.x);K=1./K*vec4(cos(T),sin(T),cos(T)*sin(T)+.5,1);}
