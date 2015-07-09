%% Q20

Q17
Q18

figure;
subplot(1, 2, 1);
surf(z); 
colormap gray
shading interp;
axis('tight');
view(110,45);
axis('off');
camlight
title('|\nablaz|=F(x,y)')
subplot(1, 2, 2);
surf(Z); 
colormap gray
shading interp;
axis('tight');
view(110,45);
axis('off');
camlight
title('|\nabla^2z|=p_x + q_y')


