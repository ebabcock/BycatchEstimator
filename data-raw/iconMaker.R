
library(hexSticker)

imgurl <- paste0(getwd(),"/man/figures/fish-g119836312_640.png")
sticker(imgurl,
        package="BycatchEstimator",
        p_size=15,
        p_color = "#212529",
        h_fill = "lightgrey",
        h_color = "grey",
        s_x=1, s_y=.75,
        s_width=.5,
        filename="man/figures/imgfile.png")
