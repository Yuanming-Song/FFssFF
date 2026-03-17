# Plot AA contact analysis: generates multi-panel time-series of
# intermolecular contacts (total, pi-pi, backbone-backbone) for
# acidic/basic solvated and vacuum FFssFF simulations.
# Sourced by plot_contact_analysis.Rmd; expects main_dir, step,
# pltthemetemp, save_plots, pltsavedir variables to be set.

# Define color scheme
color_scheme <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

# Define RDA filenames
acid_sol_rda <- paste0("FFssFF_Acid_mix_water_ions_step", step, ".rda")
base_sol_rda <- paste0("FFssFF_Base_mix_water_ions_step", step, ".rda")
base_vac_rda <- paste0("FFssFF_Base_mix_lattice_step", step, ".rda")

# Load the data
load(file.path(main_dir, "AA/Analysis/contact_analysis/output", acid_sol_rda))
acid_sol_data <- contact_data
rm(contact_data)

load(file.path(main_dir, "AA/Analysis/contact_analysis/output", base_sol_rda))
base_sol_data <- contact_data
rm(contact_data)

load(file.path(main_dir, "AA/Analysis/contact_analysis/output", base_vac_rda))
base_vac_data <- contact_data
rm(contact_data)

# Calculate common y-limits
y_limits <- list(
    range(acid_sol_data$total_contacts, acid_sol_data$pi_pi_contacts, acid_sol_data$bb_bb_contacts),
    range(base_sol_data$total_contacts, base_sol_data$pi_pi_contacts, base_sol_data$bb_bb_contacts),
    range(base_vac_data$total_contacts, base_vac_data$pi_pi_contacts, base_vac_data$bb_bb_contacts)
)

# Set common y-limits with padding
y_min <- min(sapply(y_limits, min)) * 0.95  # 5% padding at bottom
y_max <- max(sapply(y_limits, max)) * 1.05  # 5% padding at top

# Function to create contact plot
create_contact_plot <- function(data, title) {
    # Filter data by step
    data <- data[seq(1, nrow(data), step),]
    
    # Reshape data for plotting
    plot_data <- data %>%
        pivot_longer(
            cols = c(pi_pi_contacts, bb_bb_contacts, total_contacts),
            names_to = "interaction_type",
            values_to = "count"
        ) %>%
        mutate(interaction_type = factor(interaction_type,
            levels = c("total_contacts", "pi_pi_contacts", "bb_bb_contacts"),
            labels = c("Total", "Pi-Pi", "Backbone-Backbone")
        ))
    
    ggplot(plot_data, aes(x = frame, y = count, color = interaction_type)) +
        geom_line() +
        scale_color_manual(values = color_scheme[1:3], name = "Interaction Type") +
        labs(
            title = title,
            x = "Frame",
            y = "Number of Contacts"
        ) +
        ylim(y_min, y_max) +  # Set common y-limits
        pltthemetemp
}

# Create plots for each system
plots <- list()

# Acid in water
plots[["acid"]] <- create_contact_plot(
    acid_sol_data,
    "Charged FFssFF in Water"
)

# Base in water
plots[["base_sol"]] <- create_contact_plot(
    base_sol_data,
    "Neutral FFssFF in Water"
)

# Base in vacuum
plots[["base_vac"]] <- create_contact_plot(
    base_vac_data,
    "Neutral FFssFF in Vacuum"
)

# Remove individual legends and axis labels
plots <- lapply(plots, function(p) p + theme(legend.position = "none", axis.title = element_blank()))

# Get legend from one plot (positioned on top)
legend <- get_legend(plots[[1]] + theme(legend.position = "right"))

# Create combined plot
combined_plot <- plot_grid(
    plotlist = plots,
    ncol = 2,
    nrow = 2,
    align = "hv",
   # labels = LETTERS[1:length(plots)],
    label_size = 10
)

# Add common labels
combined_plot <- ggdraw() +
    draw_plot(combined_plot, x = 0.02, y = 0.02, width = 0.98, height = 0.97) +
    draw_label("Counts", x = 0.02, y = 0.5, angle = 90, size = 12) +
    draw_label("Time (ns)", x = 0.5, y = 0.03, size = 12)

# Add common legend on top
final_plot <- plot_grid(
    legend,
    combined_plot,
    nrow = 2,
    rel_heights = c(0.1, 0.9)
)

print(final_plot)

# Save the plot if save_plots is TRUE
if (save_plots) {
    ggsave(
        file.path(pltsavedir, "AA_contact_analysis_plots.pdf"),
        final_plot,
        width = plot_width,
        height = plot_height,
        dpi = plot_dpi
    )
} 
