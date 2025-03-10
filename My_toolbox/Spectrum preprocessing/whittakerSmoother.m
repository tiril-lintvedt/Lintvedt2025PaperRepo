
function z = whittakerSmoother(y,lambda) 
% Reference: "A perfect smoother", PHC Eilers - Analytical chemistry, 2003
% This is the simplest formulation - check Eilers' article for an alternative 
% implementation making use of sparse matrixes.

    m = length(y);
    I = eye(m);
    D = diff(I);
    z = (I + lambda * D' * D)\y;
end
         
% Test
% y1 = [4114. 3281. 3150. 3447. 3602. 4033. 4510. 4092. 4258. 3987. 4033. 3957. 3917. 3551. 3782. 3987. 3640. 3796. 3841. 3920. 3838. 3696. 3920.]';
% y2 = [3587. 3552. 3717. 3895. 4013. 4643. 4090. 3675. 2909. 2870. 2530. 2441. 2719. 2737. 2435. 2781. 2782. 2691. 2737. 2863. 2661. 2305. 2418.]';
% whittakerSmoother(y1, 0.1)
% whittakerSmoother(y2, 0.1)