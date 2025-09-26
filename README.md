# microglia-motility

This code takes .tif files that contain z-projected, registered, concatenated images of microglia at multiple timepoints, performs Phansalkar's local thresholding method to binarize the images, and then calculates the difference in microglia process positions between timepoints to calculate motility indices. The mean motility indices are calucalted by averaging the calculated indices between each sucessive imaging timepoint.


motility index calculated by (extendedPx + retractedPx)/(thresholdedPx)
s_motility index calculated by (extendedPx + retractedPx)/(stablePx)
extension index calculated by extendedPx/stablePx
retraction index calculated by retractedPx/stablePx
stability index calculated by (extendedPx at t_n that turned into stablePx at t_(n+1))/(extendedPx)
instability index calculated by (stablePx at t_n that turned into retractedPx at t_(n+1))/(stablePx)
