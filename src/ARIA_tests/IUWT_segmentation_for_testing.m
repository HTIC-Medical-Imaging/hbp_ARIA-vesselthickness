function bw = IUWT_segmentation_for_testing(im, levels, percent, remove, fill, dark, bw_mask)
% The IUWT segmentation algorithm.
% This is a stripped-down version of the general function implemented in
% ARIA.  It excludes some of the error checking and additional prompts, so
% is somewhat easier to read, but the steps taken are the same.
% It is intended to help get the timing of the segmentation alone for DRIVE
% database images.
%
% For a general function, use SEG_IUWT in ARIA instead.
% 
% Input:
%   IM - a 2D matrix containing the image data.
%   LEVELS - a numeric vector containing the wavelet levels that
%   should be used.
%   PERCENT - threshold defined as a proportion of the pixels in the
%   image or FOV.
%   REMOVE - the minimum size an object needs to exceed in order to
%   be kept, defined in pixels.
%   FILL - the minimum size of a 'hole' (i.e. an undetected region
%   entirely surrounded by detected pixels), defined in pixels.
%   DARK - TRUE if vessels are darker than their surroundings (e.g.
%   fundus images), FALSE if they are brighter (e.g. fluorescein
%   angiograms).
%   BW_MASK - a logical matrix, the same size as IM, representing a FOV
%   mask.  If this is not given, a mask is generated by applying a fixed
%   threshold of 20 to IM and applying a 3x3 erosion.  This is suitable for
%   the DRIVE database, but may well be inappropriate for other images.
% 
% Output:
%   BW - the binary image produced by the algorithm, in which vessel pixels
%   are TRUE and non-vessel pixels are FALSE.
% 
% 
% Note that this function is intended for use with DRIVE database images,
% and so includes optional mask generation - because this is actually one
% of the slower parts of segmentation (at least using older MATLAB
% versions) for DRIVE database images, and we did not want to over-report 
% the speed of our segmentation approach by excluding it from timing
% measurements.  However the function could also be used with other image
% sources, in which case care should be taken to set the value of BW_MASK
% appropriately.
%
%
% Copyright � 2011 Peter Bankhead.
% See the file : Copyright.m for further details.


% Extract the green channel from the image if necessary
if size(im, 3) > 1
    im = single(im(:,:,2));
end

% Create a mask, if necessary
if nargin < 7
    bw_mask = imerode(im > 20, ones(3));
end

% Apply the actual segmentation
w = iuwt_vessels(im, levels);
bw = percentage_segment(w, percent, dark, bw_mask);
bw = clean_segmented_image(bw, remove, fill);