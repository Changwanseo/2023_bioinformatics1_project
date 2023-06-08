# -*- coding: utf-8 -*-
from .GenomeTrack import GenomeTrack
import numpy as np


class PileupTrack(GenomeTrack):
    SUPPORTED_ENDINGS = [
        ".txt"
    ]  # this is used by make_tracks_file to guess the type of track based on file name
    TRACK_TYPE = "alignment"

    def __init__(
        self,
        tracks_file,
        fig_width=DEFAULT_FIGURE_WIDTH,
        fig_height=None,
        fontsize=None,
        dpi=None,
        track_label_width=0.1,
        plot_regions=None,
        plot_width=None,
    ):
        self.fig_width = fig_width
        self.fig_height = fig_height
        self.dpi = dpi
        self.type_list = None
        self.track_list = None
        start = self.print_elapsed(None)
        self.available_tracks = self.get_available_tracks()
        self.available_types = self.get_available_types()
        self.parse_tracks(tracks_file, plot_regions=plot_regions)
        if fontsize:
            fontsize = fontsize
        else:
            fontsize = float(fig_width) * 0.3
        # the track label width is the fraction of
        # the figure width that is used
        # for the track 'title' or label.
        self.width_ratios = (0.01, 1 - track_label_width, track_label_width)

        # Process the width:
        if plot_width is not None:
            self.fig_width = (
                plot_width
                / (DEFAULT_MARGINS["right"] - DEFAULT_MARGINS["left"])
                * (1 + 2 / 3 * 0.01)
                / (self.width_ratios[1] / sum(self.width_ratios))
            )

        font = {"size": fontsize}
        matplotlib.rc("font", **font)
        # initialize each track
        self.track_obj_list = []
        for idx, properties in enumerate(self.track_list):
            log.info(f"initialize {properties['section_name']}")
            track_class = properties["track_class"]
            properties["region"] = plot_regions.copy()
            self.track_obj_list.append(track_class(properties))

        # initialize each type
        self.type_obj_list = []
        for idx, properties in enumerate(self.type_list):
            log.info(f"initialize {properties['section_name']}")
            track_class = properties["track_class"]
            properties["region"] = plot_regions.copy()
            self.type_obj_list.append(track_class(properties))

        log.info("time initializing track(s):")
        self.print_elapsed(start)

    def plot(
        self,
        file_name,
        chrom,
        start,
        end,
        title=None,
        h_align_titles="left",
        decreasing_x_axis=False,
    ):
        track_height = self.get_tracks_height(start_region=start, end_region=end)

        if self.fig_height:
            fig_height = self.fig_height
        else:
            fig_height = sum(track_height) / (
                DEFAULT_MARGINS["top"] - DEFAULT_MARGINS["bottom"]
            )

        log.debug(
            f"Figure size in cm is {self.fig_width} x {fig_height}."
            f" Dpi is set to {self.dpi}\n"
        )
        fig = plt.figure(figsize=self.cm2inch(self.fig_width, fig_height))

        fig.subplots_adjust(
            wspace=0,
            hspace=0.0,
            left=DEFAULT_MARGINS["left"],
            right=DEFAULT_MARGINS["right"],
            bottom=DEFAULT_MARGINS["bottom"],
            top=DEFAULT_MARGINS["top"],
        )

        if title:
            fig.suptitle(title)

        grids = matplotlib.gridspec.GridSpec(
            len(track_height),
            3,
            height_ratios=track_height,
            width_ratios=self.width_ratios,
            wspace=0.01,
        )
        axis_list = []
        # skipped_tracks is the count of tracks that have the
        # 'overlay_previous' parameter and should be skipped
        skipped_tracks = 0
        plot_axis = None
        for idx, track in enumerate(self.track_obj_list):
            log.info(f"plotting {track.properties['section_name']}")

            if track.properties["overlay_previous"] in ["yes", "share-y"]:
                overlay = True
                skipped_tracks += 1
            else:
                overlay = False

            if track.properties["overlay_previous"] == "share-y":
                ylim = plot_axis.get_ylim()
            else:
                idx -= skipped_tracks
                plot_axis = axisartist.Subplot(fig, grids[idx, 1])
                fig.add_subplot(plot_axis)
                # turns off the lines around the tracks
                plot_axis.axis[:].set_visible(False)
                # to make the background transparent
                plot_axis.patch.set_visible(False)
                if not overlay:
                    y_axis = plt.subplot(grids[idx, 0])
                    y_axis.set_axis_off()

                    label_axis = plt.subplot(grids[idx, 2])
                    label_axis.set_axis_off()
                    # I get the width of the label_axis to be able to wrap the
                    # labels when right or center aligned.
                    width_inch = label_axis.get_window_extent().width
                    width_dpi = width_inch * self.dpi / fig.dpi

            if decreasing_x_axis:
                plot_axis.set_xlim(end, start)
            else:
                plot_axis.set_xlim(start, end)
            track.plot(plot_axis, chrom, start, end)
            track.plot_y_axis(y_axis, plot_axis)
            track.plot_label(label_axis, width_dpi=width_dpi, h_align=h_align_titles)

            if track.properties["overlay_previous"] == "share-y":
                plot_axis.set_ylim(ylim)

            if not overlay:
                axis_list.append(plot_axis)

        for current_type in self.type_obj_list:
            current_type.plot(axis_list, fig, chrom, start, end)

        fig.savefig(file_name, dpi=self.dpi, transparent=False)
        return fig
