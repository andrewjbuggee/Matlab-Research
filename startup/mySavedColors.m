% --- My Saved Colors ---

% A hodgepodge of colors I've saved over time because I found tehm
% appealing

% INPUTS
%   n - the number of colors you'd like
%   'random' - the string input tells the function that you want either
%   some precise color, or a random assortment of colors. The input
%   possibilities are:
%       (a) 'random' - function will output a random assortment of n colors
%       (b) 'fixed' - function will output the color codes associated with
%       the vector n according to how they have been saved in the matrix
%       below

function C = mySavedColors(n, randomORfixed)


savedColors = [8.492724501845559e-01     5.437108503990062e-02     9.681090252965144e-01;...        % (1) a nice bright pink
    4.843239544631898e-01     7.049687413641152e-01     6.568033076805693e-01;...                   % (2) sea foam green
    0.939001561999887   0.587044704531417   0.230488160211558;...                                   % (3) Kraft Mac and Cheese Orange
    7.065829534842050e-02     9.583877169541491e-01     6.998897395776505e-01;...                   % (4) a bright aquamarine
    0.827450980392157, 0.368627450980392, 0.376470588235294;...                                     % (5) A pleasing salmon red
    3.563645953575846e-01     4.380836048512262e-01     5.147715889386915e-01;...                   % (6) a pleasing blueish gray
    4.672877744560139e-01     8.367186588778170e-01     2.206248068838936e-01;...                   % (7) lime green
    8.259258869082566e-01     6.171862413652410e-01     1.938744022737314e-01;...                   % (8) matte harvest gold
    0.350727103576883   0.622475086001227   0.470923348517591;...                                   % (9) Matte Irish green
    0.550156342898422   0.301246330279491   0.194764289567049;...                                   % (10) UPS brown
    0.80,0.79,0.85;...                                                                              % (11) A pale grey
    0.96,0.42,0.65;...                                                                              % (12) Bubble gum pink
    0.49,0.92,0.04;...                                                                              % (13) Neon green
    0.14,0.96,0.93;                                                                                 % (14) A bright electric blue
    6.017416069441683e-02     4.714776469952955e-01     3.110165739381916e-01;...                   % (15) forest green
     ];


if strcmp(randomORfixed, 'random')==true

    indices2use = randperm(size(savedColors,1), n);          % random set of n unique integers between the values 1 and the length of my saved colors.
    C = savedColors(indices2use,:);

elseif strcmp(randomORfixed, 'fixed')==true
    
    C = savedColors(n,:);


else

    error([newline,'I dont know if you want a random assortment of colors or some precise colors', newline])

end





end