function denoised = PRI_NL_PCA(img,v,t,map,Rician)
B = padarray(img,[v, v, v],'symmetric');
s = size(B);
denoised = zeros(size(img));
for i = v+1:s(1)-v
    i-v
    for j = v+1:s(2)-v
        for l = v+1:s(3)-v
          search_blk = B(i-v:i+v, j-v:j+v, l-v:l+v);
          ref_patch = B(i-t:i+t, j-t:j+t, l-t:l+t);
          
          sig_wiegh = 0;
          sig_patch_wiegh = 0;
          h = 0.4*map(i-v,j-v,l-v);
          for k=t+1:2*v-t+1
              for q=t+1:2*v-t+1
                  for w=t+1:2*v-t+1
                      patch = search_blk(k-t:k+t,q-t:q+t,w-t:w+t);
                      
                      gi = (ref_patch(t+1,t+1, t+1)-patch(t+1,t+1, t+1))^2;
                      ui =3*(mean(ref_patch(:))-mean(patch(:)))^2;
                      
                      
                      wiegh = exp(-0.5*((gi+ui)/2*h^2));
                      sig_wiegh = sig_wiegh + wiegh;
                      if Rician
                          patch_wiegh = wiegh*(patch(t+1,t+1,t+1))^2;
                      else
                      patch_wiegh = wiegh*patch(t+1,t+1,t+1);
                      end
                      sig_patch_wiegh = sig_patch_wiegh+patch_wiegh;
                  end
              end
          end
          if Rician
              denoised(i-v,j-v,l-v) =sqrt(max( (sig_patch_wiegh/sig_wiegh)-2*(map(i-v,j-v,l-v))^2,0));
          else
          denoised(i-v,j-v,l-v) = sig_patch_wiegh/sig_wiegh;
          end
        end
    end
end
                      
end                      



