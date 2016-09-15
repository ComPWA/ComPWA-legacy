#ifndef ROOTPLOTHELPER_HPP_
#define ROOTPLOTHELPER_HPP_

#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <set>

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TLatex.h"
#include "TColor.h"
#include "TLine.h"
#include "TF1.h"

namespace NeatPlotting {

// ==================== namespace helper structures & functions ====================

static Int_t *normal_colors = 0;
static Int_t *diff_colors = 0;
static int palette_mode = -1;

inline void setDiffColorPalette(const unsigned int color_gradient_steps = 255) {
	if (diff_colors == 0) {
		Double_t Red[3] = { 0.0, 1.0, 1.0 };
		Double_t Green[3] = { 0.0, 1.0, 0.0 };
		Double_t Blue[3] = { 1.0, 1.0, 0.0 };
		Double_t Length[3] = { 0.0, 0.50, 1.0 };

		diff_colors = new Int_t[color_gradient_steps];
		Int_t FI = TColor::CreateGradientColorTable(3, Length, Red, Green, Blue,
				color_gradient_steps);
		for (unsigned int i = 0; i < color_gradient_steps; i++)
			diff_colors[i] = FI + i;
	}

	if (1 != palette_mode) {
		gROOT->UseCurrentStyle();
		gStyle->SetPalette(color_gradient_steps, diff_colors);
		gStyle->SetNumberContours(color_gradient_steps);
		palette_mode = 1;
	}
}

inline void setNormalColorPalette(
		const unsigned int color_gradient_steps = 255) {
	const Int_t NRGBs = 5;

	if (normal_colors == 0) {

		Double_t Red[NRGBs] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
		Double_t Green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
		Double_t Blue[NRGBs] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
		Double_t Length[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };

		normal_colors = new Int_t[color_gradient_steps];
		Int_t FI = TColor::CreateGradientColorTable(NRGBs, Length, Red, Green, Blue,
				color_gradient_steps);
		for (unsigned int i = 0; i < color_gradient_steps; i++)
			normal_colors[i] = FI + i;
	}
	if (0 != palette_mode) {
		gROOT->UseCurrentStyle();
		gStyle->SetPalette(color_gradient_steps, normal_colors);
		gStyle->SetNumberContours(color_gradient_steps);
		palette_mode = 0;
	}
}

template<class T> class DrawableDataObjectDrawOptionPair {
	friend class PlotBundle;

	std::string draw_option;

	void removeAxisDrawOption() {
		char chars[] = "aA";

		for (unsigned int i = 0; i < strlen(chars); ++i) {
			// you need include <algorithm> to use general algorithms like std::remove()
			draw_option.erase(
					std::remove(draw_option.begin(), draw_option.end(), chars[i]),
					draw_option.end());
		}
	}

public:
	T data_object;

	DrawableDataObjectDrawOptionPair(T data_object_, std::string draw_option_) :
			data_object(data_object_), draw_option(draw_option_) {
		removeAxisDrawOption();
	}
};

// style objects

struct MarkerStyle {
	int marker_style;
	int marker_color;
	double marker_size;

	MarkerStyle() :
			marker_style(24), marker_color(1), marker_size(0.6) {
	}
};

struct LineStyle {
	int line_style;
	int line_color;
	int line_width;

	LineStyle() :
			line_style(1), line_color(1), line_width(1) {
	}
};

struct DataObjectStyle {
	MarkerStyle marker_style;
	LineStyle line_style;
	std::string draw_option;
};

struct TextStyle {
	int text_font;
	double text_size;
	int text_color;

	TextStyle() :
			text_font(132), text_size(0.05), text_color(1) {
	}

	bool operator<(const TextStyle &rhs) const {
		if (text_font < rhs.text_font)
			return true;
		else if (text_font > rhs.text_font)
			return false;
		if (text_size < rhs.text_size)
			return true;
		else if (text_size > rhs.text_size)
			return false;
		if (text_color < rhs.text_color)
			return true;
		else if (text_color > rhs.text_color)
			return false;
		return false;
	}
	bool operator>(const TextStyle &rhs) const {
		return (rhs < *this);
	}
};

struct AxisStyle {
	TextStyle axis_text_style;

	bool log_scale;

	// axis styles
	double axis_label_text_offset;
	double axis_title_text_offset;

	AxisStyle() :
			log_scale(false), axis_label_text_offset(0.005), axis_title_text_offset(
					1.15) {
	}

	bool operator<(const AxisStyle &rhs) const {
		if (axis_text_style < rhs.axis_text_style)
			return true;
		else if (axis_text_style > rhs.axis_text_style)
			return false;
		if (log_scale < rhs.log_scale)
			return true;
		else if (log_scale > rhs.log_scale)
			return false;
		if (axis_label_text_offset < rhs.axis_label_text_offset)
			return true;
		else if (axis_label_text_offset > rhs.axis_label_text_offset)
			return false;
		if (axis_title_text_offset < rhs.axis_title_text_offset)
			return true;
		else if (axis_title_text_offset > rhs.axis_title_text_offset)
			return false;
		return false;
	}
	bool operator>(const AxisStyle &rhs) const {
		return (rhs < *this);
	}
};

struct PlotStyle {
	AxisStyle x_axis_style;
	AxisStyle y_axis_style;
	AxisStyle z_axis_style;

	// other
	int statistics_options; // 0 none higher see setoptstat in root
	int palette_color_style; // 0 normal, 1 diff/ratio

	PlotStyle() :
			palette_color_style(0), statistics_options(0) {
	}

	bool operator<(const PlotStyle &rhs) const {
		if (x_axis_style < rhs.x_axis_style)
			return true;
		else if (x_axis_style > rhs.x_axis_style)
			return false;
		if (y_axis_style < rhs.y_axis_style)
			return true;
		else if (y_axis_style > rhs.y_axis_style)
			return false;
		if (z_axis_style < rhs.z_axis_style)
			return true;
		else if (z_axis_style > rhs.z_axis_style)
			return false;
		if (statistics_options < rhs.statistics_options)
			return true;
		else if (statistics_options > rhs.statistics_options)
			return false;
		if (palette_color_style < rhs.palette_color_style)
			return true;
		else if (palette_color_style > rhs.palette_color_style)
			return false;
		return false;
	}
	bool operator>(const PlotStyle &rhs) const {
		return (rhs < *this);
	}

};

// data structures and wrappers

struct GraphPoint {
	double x;
	double y;
	double x_err_low;
	double x_err_high;
	double y_err_low;
	double y_err_high;
	GraphPoint() :
			x(0.0), y(0.0), x_err_low(0.0), x_err_high(0.0), y_err_low(0.0), y_err_high(
					0.0) {
	}
};

struct PadBoundaries {
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;

	bool is_2d;
	bool is_empty;

	PadBoundaries() :
			x_min(0.0), x_max(0.0), y_min(0.0), y_max(0.0), is_empty(true), is_2d(
					false) {
	}

	PadBoundaries(TH1 &hist) :
			is_empty(false), is_2d(false) {
		x_min = hist.GetXaxis()->GetXmin();
		x_max = hist.GetXaxis()->GetXmax();

		if (hist.GetDimension() == 1) {
			if (gPad) {
				gPad->Update();
				if (gPad->GetLogy()) {
					y_min = pow(10, gPad->GetUymin());
					y_max = pow(10, gPad->GetUymax());
				} else {
					y_min = gPad->GetUymin();
					y_max = gPad->GetUymax();
				}
			}
		}

		if (hist.GetDimension() > 1) {
			is_2d = true;
			y_min = hist.GetYaxis()->GetXmin();
			y_max = hist.GetYaxis()->GetXmax();
			z_min = hist.GetMinimum();
			z_max = hist.GetMaximum();
		}
	}
	PadBoundaries getUnionBoundaries(PadBoundaries & other_boundaries) {
		PadBoundaries union_pad;
		if (!is_empty && !other_boundaries.is_empty) {
			union_pad.is_empty = false;
			if (x_min < other_boundaries.x_min)
				union_pad.x_min = x_min;
			else
				union_pad.x_min = other_boundaries.x_min;
			if (x_max > other_boundaries.x_max)
				union_pad.x_max = x_max;
			else
				union_pad.x_max = other_boundaries.x_max;
			if (y_min < other_boundaries.y_min)
				union_pad.y_min = y_min;
			else
				union_pad.y_min = other_boundaries.y_min;
			if (y_max > other_boundaries.y_max)
				union_pad.y_max = y_max;
			else
				union_pad.y_max = other_boundaries.y_max;
		}

		return union_pad;
	}
};

struct AxisRange {
	double low;
	double high;
	bool active;

	AxisRange() :
			active(false) {
	}
	AxisRange(double low_, double high_) :
			low(low_), high(high_), active(true) {
	}
};

inline double calculateAbsolutePositionFromPadBottom(
		double relative_position_from_bottom) {
	double abs_pos_from_bottom(0.0);

	if (gPad->GetLogy()) {
		abs_pos_from_bottom = pow(10,
				gPad->GetUymin()
						+ (gPad->GetUymax() - gPad->GetUymin())
								* relative_position_from_bottom);
	} else {
		abs_pos_from_bottom =
		gPad->GetUymin()
				+ (gPad->GetUymax() - gPad->GetUymin()) * relative_position_from_bottom;
	}
	return abs_pos_from_bottom;
}

inline double calculateAbsolutePositionFromPadLeft(
		double relative_position_from_left) {
	double abs_pos_from_left(0.0);

	if (gPad->GetLogx()) {
		abs_pos_from_left = pow(10,
				gPad->GetUxmin()
						+ (gPad->GetUxmax() - gPad->GetUxmin())
								* relative_position_from_left);
	} else {
		abs_pos_from_left =
		gPad->GetUxmin()
				+ (gPad->GetUxmax() - gPad->GetUxmin()) * relative_position_from_left;
	}
	return abs_pos_from_left;
}

class PlotLabel {
	friend class PlotBundle;
	TLatex label;
	double relative_position_from_bottom;
	double relative_position_from_left;
	double absolute_position_x;
	double absolute_position_y;
	bool use_relative_positioning;

public:
	PlotLabel(const std::string &text_, const TextStyle &text_style) :
			label(0.0, 0.0, text_.c_str()), relative_position_from_bottom(0.0), relative_position_from_left(
					0.0), absolute_position_x(0.0), absolute_position_y(0.0), use_relative_positioning(
					true) {
		label.SetTextColor(text_style.text_color);
		label.SetTextSize(text_style.text_size);
		label.SetTextFont(text_style.text_font);
	}
	PlotLabel(const PlotLabel &rhs) :
			label(0.0, 0.0, rhs.label.GetTitle()), relative_position_from_bottom(
					rhs.relative_position_from_bottom), relative_position_from_left(
					rhs.relative_position_from_left), absolute_position_x(
					rhs.absolute_position_x), absolute_position_y(
					rhs.absolute_position_y), use_relative_positioning(
					rhs.use_relative_positioning) {
		label.SetTextColor(rhs.label.GetTextColor());
		label.SetTextSize(rhs.label.GetTextSize());
		label.SetTextFont(rhs.label.GetTextFont());
	}

	PlotLabel& operator=(const PlotLabel &rhs) {
		if (this != &rhs) {
			label.SetTitle(rhs.label.GetTitle());
			label.SetTextColor(rhs.label.GetTextColor());
			label.SetTextSize(rhs.label.GetTextSize());
			label.SetTextFont(rhs.label.GetTextFont());

			absolute_position_x = rhs.absolute_position_x;
			absolute_position_y = rhs.absolute_position_y;
			relative_position_from_bottom = rhs.relative_position_from_bottom;
			relative_position_from_left = rhs.relative_position_from_left;
		}
		return *this;
	}

	std::string getTitle() {
		return label.GetTitle();
	}

	void setTitle(std::string title_) {
		label.SetTitle(title_.c_str());
	}

	void setAbsolutionPosition(double x, double y) {
		absolute_position_x = x;
		absolute_position_y = y;
		use_relative_positioning = false;
	}
	void setRelativePosition(double x, double y) {
		relative_position_from_bottom = y;
		relative_position_from_left = x;
		absolute_position_x = true;
	}

	void calculateAbsolutePositioning() {
		if (use_relative_positioning) {
			label.SetX(
					calculateAbsolutePositionFromPadLeft(relative_position_from_left));
			label.SetY(
					calculateAbsolutePositionFromPadBottom(
							relative_position_from_bottom));
		} else {
			label.SetX(absolute_position_x);
			label.SetY(absolute_position_y);
		}
	}
};

struct PlotAccessories {
	bool draw_labels;

	// todo: this should go into some layout manager class
	bool use_line_layout;
	double label_text_leftpos;
	double label_text_toppos;
	double label_text_spacing;

	std::vector<PlotLabel> labels;

	std::vector<double> x_parallel_lines;
	std::vector<double> y_parallel_lines;

	PlotAccessories() :
			draw_labels(true), use_line_layout(false) {
	}
};

class PlotAxis {
	TH2D* axis_dummy;

	void beautifyAxis(TAxis *axis, const AxisStyle &style) const {
		// axis settings
		axis->SetTitleFont(style.axis_text_style.text_font);
		axis->SetLabelFont(style.axis_text_style.text_font);
		axis->SetLabelSize(style.axis_text_style.text_size);
		axis->SetTitleSize(style.axis_text_style.text_size);
		axis->SetLabelOffset(style.axis_label_text_offset);
		axis->SetTitleOffset(style.axis_title_text_offset);
	}

public:
	std::string x_axis_title;
	std::string y_axis_title;
	std::string z_axis_title;

	AxisRange x_axis_range;
	AxisRange y_axis_range;
	AxisRange z_axis_range;

	PlotAxis() :
			axis_dummy(0) {
	}
	virtual ~PlotAxis() {
		/*if (axis_dummy) {
		 std::cout<<"aasdf1"<<std::endl;
		 delete (axis_dummy);
		 std::cout<<"aasdf2"<<std::endl;
		 }*/
	}

	void initAxis(const PadBoundaries &pad_boundaries) {
		if (!x_axis_range.active) {
			x_axis_range.low = pad_boundaries.x_min;
			x_axis_range.high = pad_boundaries.x_max;
		}

		if (!y_axis_range.active) {
			y_axis_range.low = pad_boundaries.y_min;
			y_axis_range.high = pad_boundaries.y_max;
		}

		if (!z_axis_range.active) {
			z_axis_range.low = pad_boundaries.z_min;
			z_axis_range.high = pad_boundaries.z_max;
		}

		axis_dummy = new TH2D("axis_hist", "", 1, x_axis_range.low,
				x_axis_range.high, 1, y_axis_range.low, y_axis_range.high);
		axis_dummy->Fill(0.0, 0.0, 0.0);

		if (pad_boundaries.is_2d) {
			axis_dummy->GetZaxis()->SetLimits(z_axis_range.low, z_axis_range.high);
			axis_dummy->GetZaxis()->SetRangeUser(z_axis_range.low, z_axis_range.high);

			axis_dummy->Draw("COLZ");
		} else
			axis_dummy->Draw("");

		axis_dummy->GetXaxis()->SetTitle(x_axis_title.c_str());
		axis_dummy->GetYaxis()->SetTitle(y_axis_title.c_str());
		axis_dummy->GetZaxis()->SetTitle(z_axis_title.c_str());
	}

	void redrawAxis() const {
		axis_dummy->Draw("SAMEAXIS");
		gPad->Update();
		gPad->RedrawAxis();
	}

	void beautifyAxis(const PlotStyle &style) {
		axis_dummy->SetStats(style.statistics_options);
		axis_dummy->SetTitle("");

		beautifyAxis(axis_dummy->GetXaxis(), style.x_axis_style);
		beautifyAxis(axis_dummy->GetYaxis(), style.y_axis_style);
		beautifyAxis(axis_dummy->GetZaxis(), style.z_axis_style);

		// set logarithmic scale
		if (style.x_axis_style.log_scale) {
			if (axis_dummy->GetXaxis()->GetXmin() >= 0.0) {
				if (axis_dummy->GetXaxis()->GetXmin() == 0.0) {
					axis_dummy->GetXaxis()->SetRangeUser(1e-1,
							axis_dummy->GetXaxis()->GetXmax());
					axis_dummy->GetXaxis()->SetLimits(1e-1,
							axis_dummy->GetXaxis()->GetXmax());
				}
				gPad->SetLogx(1);
			} else
				std::cout << "log scale for negative axis range not allowed!"
						<< std::endl;
		} else
			gPad->SetLogx(0);
		if (style.y_axis_style.log_scale) {
			if (axis_dummy->GetYaxis()->GetXmin() >= 0.0) {
				if (axis_dummy->GetYaxis()->GetXmin() == 0.0) {
					axis_dummy->GetYaxis()->SetRangeUser(1e-1,
							axis_dummy->GetYaxis()->GetXmax());
					axis_dummy->GetYaxis()->SetLimits(1e-1,
							axis_dummy->GetYaxis()->GetXmax());
				}
				gPad->SetLogy(1);
			} else
				std::cout << "log scale for negative axis range not allowed!"
						<< std::endl;
		} else
			gPad->SetLogy(0);
		if (style.z_axis_style.log_scale) {
			if (axis_dummy->GetZaxis()->GetXmin() >= 0.0) {
				if (axis_dummy->GetZaxis()->GetXmin() == 0.0) {
					axis_dummy->GetZaxis()->SetRangeUser(1e-1,
							axis_dummy->GetZaxis()->GetXmax());
					axis_dummy->GetZaxis()->SetLimits(1e-1,
							axis_dummy->GetZaxis()->GetXmax());
				}
				gPad->SetLogz(1);
			} else
				std::cout << "log scale for negative axis range not allowed!"
						<< std::endl;
		} else
			gPad->SetLogz(0);

		if (1 == style.palette_color_style) {
			setDiffColorPalette();

			double min = axis_dummy->GetMinimum();
			double max = axis_dummy->GetMaximum();

			double maxval = std::max(-min, max);
			axis_dummy->GetZaxis()->SetLimits(-maxval, maxval);
			axis_dummy->GetZaxis()->SetRangeUser(-maxval, maxval);
		} else {
			setNormalColorPalette();
		}
	}
};

class PlotBundle {
	std::vector<DrawableDataObjectDrawOptionPair<TH1*> > histograms;
	std::vector<DrawableDataObjectDrawOptionPair<TGraph*> > graphs;
	std::vector<DrawableDataObjectDrawOptionPair<TGraph2D*> > graphs2d;

	bool drawHistsAndGraphs(const PlotStyle &plot_style) const {
		if (hasDrawables()) {
			std::vector<DrawableDataObjectDrawOptionPair<TH1*> >::const_iterator hist_it;

			for (hist_it = histograms.begin(); hist_it != histograms.end();
					hist_it++) {
				TH1* hist = hist_it->data_object;
				TString draw_option = hist_it->draw_option;
				if (hist) {
					if (hist->GetDimension() > 1) {
						hist->GetZaxis()->SetLimits(plot_axis.z_axis_range.low,
								plot_axis.z_axis_range.high);
						hist->GetZaxis()->SetRangeUser(plot_axis.z_axis_range.low,
								plot_axis.z_axis_range.high);
					}
					draw_option += "SAME";
					hist->Draw(draw_option);
				}
			}

			std::vector<DrawableDataObjectDrawOptionPair<TGraph*> >::const_iterator graph_it;

			for (graph_it = graphs.begin(); graph_it != graphs.end(); graph_it++) {
				TGraph* graph = graph_it->data_object;
				TString draw_option = graph_it->draw_option;
				if (graph) {
					draw_option += "SAME";
					graph->Draw(draw_option);
				}
			}

			std::vector<DrawableDataObjectDrawOptionPair<TGraph2D*> >::const_iterator graph2d_it;

			for (graph2d_it = graphs2d.begin(); graph2d_it != graphs2d.end();
					graph2d_it++) {
				TGraph2D* graph = graph2d_it->data_object;
				TString draw_option = graph2d_it->draw_option;
				if (graph) {
					draw_option += "SAME";
					graph->Draw(draw_option);
				}
			}
			return true;
		} else {
			std::cout
					<< "Dude, you forgot to add a histogram or graph in the plot bundle..."
					<< " Please add at least one and make sure it points to an existing object!"
					<< std::endl;
			return false;
		}
	}

public:
	PlotAccessories plot_decoration;
	PlotAxis plot_axis;
	PadBoundaries total_boundaries;

	PlotBundle() {
	}
	virtual ~PlotBundle() {
	}

	bool hasDrawables() const {
		return !total_boundaries.is_empty;
	}

	std::vector<DrawableDataObjectDrawOptionPair<TH1*> > getHistograms() const {
		return histograms;
	}
	std::vector<DrawableDataObjectDrawOptionPair<TGraph*> > getGraphs() const {
		return graphs;
	}
	std::vector<DrawableDataObjectDrawOptionPair<TGraph2D*> > getGraphs2d() const {
		return graphs2d;
	}

	void addHistogram(TH1* hist, const DataObjectStyle &style) {
		applyMarkerStyleToDrawableObject(hist, style.marker_style);
		applyLineStyleToDrawableObject(hist, style.line_style);
		histograms.push_back(
				DrawableDataObjectDrawOptionPair<TH1*>(hist, style.draw_option));

	}
	void addGraph(TGraph* graph, const DataObjectStyle &style) {
		applyMarkerStyleToDrawableObject(graph, style.marker_style);
		applyLineStyleToDrawableObject(graph, style.line_style);
		graphs.push_back(
				DrawableDataObjectDrawOptionPair<TGraph*>(graph, style.draw_option));
	}
	void addGraph(TGraph2D* graph, const DataObjectStyle &style) {
		applyMarkerStyleToDrawableObject(graph, style.marker_style);
		applyLineStyleToDrawableObject(graph, style.line_style);
		graphs2d.push_back(
				DrawableDataObjectDrawOptionPair<TGraph2D*>(graph, style.draw_option));
	}

	PadBoundaries calculatePadBoundaries() const {
		PadBoundaries union_pad_boundaries;
		bool nothing_drawn(true);

		std::vector<DrawableDataObjectDrawOptionPair<TH1*> >::const_iterator hist_it;

		for (hist_it = histograms.begin(); hist_it != histograms.end(); hist_it++) {
			TH1* hist = hist_it->data_object;
			TString draw_option(hist_it->draw_option);
			if (hist) {
				hist->Draw(draw_option);
				if (plot_axis.x_axis_range.active)
					hist->GetXaxis()->SetRangeUser(plot_axis.x_axis_range.low,
							plot_axis.x_axis_range.high);
				PadBoundaries current_pad_boundaries(*hist);

				if (nothing_drawn) {
					union_pad_boundaries = current_pad_boundaries;
					nothing_drawn = false;
				} else {
					union_pad_boundaries = union_pad_boundaries.getUnionBoundaries(
							current_pad_boundaries);
				}
			}
		}

		std::vector<DrawableDataObjectDrawOptionPair<TGraph*> >::const_iterator graph_it;

		for (graph_it = graphs.begin(); graph_it != graphs.end(); graph_it++) {
			TGraph* graph = graph_it->data_object;
			TString draw_option("A");
			draw_option += graph_it->draw_option;
			if (graph) {
				graph->Draw(draw_option);
				if (plot_axis.x_axis_range.active)
					graph->GetXaxis()->SetRangeUser(plot_axis.x_axis_range.low,
							plot_axis.x_axis_range.high);
				PadBoundaries current_pad_boundaries(*graph->GetHistogram());
				if (nothing_drawn) {
					union_pad_boundaries = current_pad_boundaries;
					nothing_drawn = false;
				} else {
					union_pad_boundaries = union_pad_boundaries.getUnionBoundaries(
							current_pad_boundaries);
				}
			}
		}

		std::vector<DrawableDataObjectDrawOptionPair<TGraph2D*> >::const_iterator graph2d_it;

		for (graph2d_it = graphs2d.begin(); graph2d_it != graphs2d.end();
				graph2d_it++) {
			TGraph2D* graph = graph2d_it->data_object;
			TString draw_option("A");
			draw_option += graph2d_it->draw_option;
			if (graph) {
				graph->Draw(draw_option);
				PadBoundaries current_pad_boundaries(*graph->GetHistogram());
				current_pad_boundaries.is_2d = true;
				current_pad_boundaries.z_min = graph->GetZmin();
				current_pad_boundaries.z_max = graph->GetZmax();
				if (nothing_drawn) {
					union_pad_boundaries = current_pad_boundaries;
					nothing_drawn = false;
				} else {
					union_pad_boundaries = union_pad_boundaries.getUnionBoundaries(
							current_pad_boundaries);
				}
			}
		}

		return union_pad_boundaries;
	}

	void initAxis(const PlotStyle &plot_style) {
		total_boundaries = calculatePadBoundaries();
		if (hasDrawables()) {
			plot_axis.initAxis(total_boundaries);
			plot_axis.beautifyAxis(plot_style);
		}
	}

	void applyMarkerStyleToDrawableObject(TAttMarker *data_object,
			const MarkerStyle &style) const {
		data_object->SetMarkerStyle(style.marker_style);
		data_object->SetMarkerColor(style.marker_color);
		data_object->SetMarkerSize(style.marker_size);
	}

	void applyLineStyleToDrawableObject(TAttLine *data_object,
			const LineStyle &style) const {
		data_object->SetLineStyle(style.line_style);
		data_object->SetLineColor(style.line_color);
		data_object->SetLineWidth(style.line_width);
	}

	void drawPlotAccessories() {
		gPad->Update();

		// draw parallel lines
		for (unsigned int i = 0; i < plot_decoration.x_parallel_lines.size(); i++) {
			TLine *line = new TLine(gPad->GetUxmin(),
					plot_decoration.x_parallel_lines[i],
					gPad->GetUxmax(), plot_decoration.x_parallel_lines[i]);
			line->Draw();
		}
		for (unsigned int i = 0; i < plot_decoration.y_parallel_lines.size(); i++) {
			TLine *line = new TLine(plot_decoration.y_parallel_lines[i],
			gPad->GetUymin(), plot_decoration.y_parallel_lines[i],
			gPad->GetUymax());
			line->Draw();
		}

		// draw labels
		if (plot_decoration.draw_labels) {
			for (unsigned int i = 0; i < plot_decoration.labels.size(); i++) {
				PlotLabel &plot_label = plot_decoration.labels[i];
				if (plot_decoration.use_line_layout) {
					// set relative position
					plot_label.setRelativePosition(plot_decoration.label_text_leftpos,
							plot_decoration.label_text_toppos
									- plot_decoration.label_text_spacing * i);
				}
				plot_label.calculateAbsolutePositioning();
				plot_label.label.Draw();
			}
		}
	}

	void drawOnCurrentPad(const PlotStyle &plot_style) {
		if (gPad) {
			initAxis(plot_style);

			bool something_was_drawn = drawHistsAndGraphs(plot_style);

			// if nothing was drawn just skip the decoration part
			if (something_was_drawn) {
				drawPlotAccessories();
			}
			plot_axis.redrawAxis();
		}
	}
};

typedef std::pair<PlotBundle, const PlotStyle*> CompletePlot;
typedef std::pair<int, int> SubpadCoordinates;
typedef std::map<SubpadCoordinates, CompletePlot> BookyPage;

class Booky {
	// all used plot bundles and styles are managed here
	std::set<PlotStyle> plot_styles;
	//std::vector<PlotBundle> plot_bundles;

	std::vector<BookyPage> booky_data;
	BookyPage current_booky_page;

	SubpadCoordinates getLargestSubpadCoordinates(
			const BookyPage &booky_page) const {
		SubpadCoordinates max_pad_coord(1, 1);

		for (std::map<SubpadCoordinates, CompletePlot>::const_iterator it =
				booky_page.begin(); it != booky_page.end(); it++) {
			SubpadCoordinates pad_coord = it->first;
			if (pad_coord.first > max_pad_coord.first)
				max_pad_coord.first = pad_coord.first;
			if (pad_coord.second > max_pad_coord.second)
				max_pad_coord.second = pad_coord.second;
		}
		return max_pad_coord;
	}

	void fillBookyPage(BookyPage &booky_page) const {
// assume the current pad is the top canvas that is divided appropriately
		TVirtualPad *current_canvas = gPad;

		gPad->Clear();

		SubpadCoordinates max_pad_coord = getLargestSubpadCoordinates(booky_page);
		current_canvas->Divide(max_pad_coord.first, max_pad_coord.second);

		for (std::map<SubpadCoordinates, CompletePlot>::iterator it =
				booky_page.begin(); it != booky_page.end(); it++) {
			SubpadCoordinates pad_coord = it->first;

			current_canvas->cd(
					max_pad_coord.first * (pad_coord.second - 1) + pad_coord.first);
			it->second.first.drawOnCurrentPad(*it->second.second);
		}
		current_canvas->cd();
	}

public:
	Booky() {
	}

	void addCurrentPageToBooky() {
		booky_data.push_back(current_booky_page);
		current_booky_page.clear();
	}

	void addPlotToCurrentBookyPage(const PlotBundle &plot_bundle,
			const PlotStyle &plot_style, SubpadCoordinates canvas_position) {
		// first insert plot style into set if not already existent
		const PlotStyle &managed_plot_style =
				*(plot_styles.insert(plot_style).first);
		current_booky_page[canvas_position] = std::make_pair(plot_bundle,
				&managed_plot_style);
	}

	void createBooky(const TString filename) {
		std::cout << "attempting to fill canvas with " << booky_data.size()
				<< " pages (filename: " << filename << ")." << std::endl;

		TCanvas c("booky", "booky", 1000, 700);

		c.Print(filename + "["); // No actual print, just open file
		for (unsigned int i = 0; i < booky_data.size(); i++) {
			c.cd(0);
			fillBookyPage(booky_data[i]);
			c.Print(filename); // actually print canvas to file
		}
		c.Print(filename + "]");
	}
};

class GraphAndHistogramHelper {
public:
	GraphAndHistogramHelper() {
	}

// ==================== graph helper functions ====================

	TGraphAsymmErrors* makeGraph(const std::vector<GraphPoint> &data) const {
		TGraphAsymmErrors *graph = new TGraphAsymmErrors(data.size());

		for (unsigned int i = 0; i < data.size(); i++) {
			const GraphPoint &gp = data[i];
			graph->SetPoint(i, gp.x, gp.y);
			graph->SetPointError(i, gp.x_err_low, gp.x_err_high, gp.y_err_low,
					gp.y_err_high);
		}

		graph->SetTitle("");
		return graph;
	}

	TGraphAsymmErrors* makeDifferenceGraph(const TGraphAsymmErrors *g1,
			const TGraphAsymmErrors *g2) const {
		if (g1->GetN() != g2->GetN())
			std::cout << "Warning: the two graphs have a different number of events!"
					<< std::endl;

		std::vector<GraphPoint> points;
		points.reserve(g1->GetN());

		for (unsigned int i = 0; i < g1->GetN(); i++) {
			double x1;
			double y1;
			g1->GetPoint(i, x1, y1);
			double x1errup = g1->GetErrorXhigh(i);
			double x1errlow = g1->GetErrorXlow(i);
			double y1errup = g1->GetErrorYhigh(i);
			double y1errlow = g1->GetErrorYlow(i);

			for (unsigned int j = 0; j < g2->GetN(); j++) {
				double x2;
				double y2;
				g2->GetPoint(j, x2, y2);

				if (x1 == x2) {
					double x2errup = g2->GetErrorXhigh(i);
					double x2errlow = g2->GetErrorXlow(i);
					double y2errup = g2->GetErrorYhigh(i);
					double y2errlow = g2->GetErrorYlow(i);

					GraphPoint gp;
					gp.x = x1;
					gp.y = y1 - y2;
					gp.x_err_low = sqrt(pow(x1errlow, 2.0) + pow(x2errlow, 2.0));
					gp.x_err_high = sqrt(pow(x1errlow, 2.0) + pow(x2errlow, 2.0));
					gp.y_err_low = sqrt(pow(y1errlow, 2.0) + pow(y2errlow, 2.0));
					gp.y_err_high = sqrt(pow(y1errlow, 2.0) + pow(y2errlow, 2.0));
					points.push_back(gp);
					break;
				}
			}
		}
		if (points.capacity() > points.size())
			points.swap(points);
		return makeGraph(points);
	}

	TGraphAsymmErrors* makeResidualZeroGraph(const TGraphAsymmErrors *g1) const {
		std::vector<GraphPoint> points;
		points.reserve(g1->GetN());

		for (unsigned int i = 0; i < g1->GetN(); i++) {
			double x1;
			double y1;
			g1->GetPoint(i, x1, y1);
			double x1errup = g1->GetErrorXhigh(i);
			double x1errlow = g1->GetErrorXlow(i);
			double y1errup = g1->GetErrorYhigh(i);
			double y1errlow = g1->GetErrorYlow(i);

			GraphPoint gp;
			gp.x = x1;
			gp.y = 0.0;
			gp.x_err_low = x1errlow;
			gp.x_err_high = x1errlow;
			gp.y_err_low = y1errlow;
			gp.y_err_high = y1errlow;
			points.push_back(gp);
		}

		return makeGraph(points);
	}

	TGraphErrors* makeDifferenceGraph(const TGraphErrors *g1,
			const TGraphErrors *g2) const {
		if (g1->GetN() != g2->GetN())
			std::cout << "Warning: the two graphs have a different number of events!"
					<< std::endl;

		std::vector<GraphPoint> points;
		points.reserve(g1->GetN());

		for (unsigned int i = 0; i < g1->GetN(); i++) {
			double x1;
			double y1;
			g1->GetPoint(i, x1, y1);
			double x1err = g1->GetErrorX(i);
			double y1err = g1->GetErrorY(i);

			for (unsigned int j = 0; j < g2->GetN(); j++) {
				double x2;
				double y2;
				g2->GetPoint(j, x2, y2);

				if (x1 == x2) {
					double x2err = g2->GetErrorX(i);
					double y2err = g2->GetErrorY(i);

					GraphPoint gp;
					gp.x = x1;
					gp.y = y1 - y2;
					gp.x_err_low = sqrt(pow(x1err, 2.0) + pow(x2err, 2.0));
					gp.x_err_high = gp.x_err_low;
					gp.y_err_low = sqrt(pow(y1err, 2.0) + pow(y2err, 2.0));
					gp.y_err_high = y1err + y2err;
					points.push_back(gp);
					break;
				}
			}
		}
		if (points.capacity() > points.size())
			points.swap(points);
		return (TGraphErrors*) makeGraph(points);
	}

	TGraphAsymmErrors* makeDifferenceGraph(const TH1D *data,
			const TGraphAsymmErrors* model) const {
		std::vector<GraphPoint> points;
// to speed things up reserve memory
		points.reserve(data->GetNbinsX());

		for (unsigned int i = 0; i < data->GetNbinsX(); i++) {
			GraphPoint gp;

			gp.x = data->GetBinCenter(i);
			gp.y = data->GetBinContent(i) - model->Eval(gp.x);
			gp.x_err_low = 0.5 * data->GetBinWidth(i);
			gp.x_err_high = 0.5 * data->GetBinWidth(i);
			gp.y_err_low = data->GetBinError(i);
			gp.y_err_high = data->GetBinError(i);

			points.push_back(gp);
		}

		return makeGraph(points);
	}

	TGraphAsymmErrors* makeRatioGraph(const TGraphAsymmErrors *g1,
			const TGraphAsymmErrors *g2) const {
		if (g1->GetN() != g2->GetN())
			std::cout << "Warning: the two graphs have a different number of points!"
					<< std::endl;

		std::vector<GraphPoint> points;
		points.reserve(g1->GetN());

		for (unsigned int i = 0; i < g1->GetN(); i++) {
			double x1;
			double y1;
			g1->GetPoint(i, x1, y1);
			double x1errup = g1->GetErrorXhigh(i);
			double x1errlow = g1->GetErrorXlow(i);
			double y1errup = g1->GetErrorYhigh(i);
			double y1errlow = g1->GetErrorYlow(i);

			for (unsigned int j = 0; j < g2->GetN(); j++) {
				double x2;
				double y2;
				g2->GetPoint(j, x2, y2);

				if (x1 == x2) {
					double x2errup = g2->GetErrorXhigh(i);
					double x2errlow = g2->GetErrorXlow(i);
					double y2errup = g2->GetErrorYhigh(i);
					double y2errlow = g2->GetErrorYlow(i);

					if (y2 != 0.0) {
						GraphPoint gp;
						gp.x = x1;
						gp.y = y1 / y2;
						gp.x_err_low = (x1errlow + x2errlow) / 2;
						gp.x_err_high = (x1errup + x2errup) / 2;
						gp.y_err_low = sqrt(
								pow(y1errlow / y2, 2.0) + pow(y2errlow * y1 / y2 / y2, 2.0));
						gp.y_err_high = sqrt(
								pow(y1errup / y2, 2.0) + pow(y2errup * y1 / y2 / y2, 2.0));
						points.push_back(gp);
					}
					break;
				}
			}
		}
		if (points.capacity() > points.size())
			points.swap(points);
		return makeGraph(points);
	}

	TGraphAsymmErrors* rescaleAxis(TGraphAsymmErrors* graph,
			double scale_factor_x) {
		std::vector<GraphPoint> points;
		points.reserve(graph->GetN());

		for (unsigned int i = 0; i < graph->GetN(); i++) {
			double x1;
			double y1;
			graph->GetPoint(i, x1, y1);
			double x1errup = graph->GetErrorXhigh(i);
			double x1errlow = graph->GetErrorXlow(i);
			double y1errup = graph->GetErrorYhigh(i);
			double y1errlow = graph->GetErrorYlow(i);

			GraphPoint gp;
			gp.x = scale_factor_x * x1;
			gp.y = y1;
			gp.x_err_low = scale_factor_x * x1errlow;
			gp.x_err_high = scale_factor_x * x1errup;
			gp.y_err_low = y1errlow;
			gp.y_err_high = y1errup;
			points.push_back(gp);
		}
		return makeGraph(points);
	}

// ==================== histogram helper functions ====================

	TH2D* makeDifferenceHistogram(TH2D *h1, TH2D *h2) const {
		TH2D* diff(0);
		if (h1 && h2) {
			diff = new TH2D(*h1);
			diff->Add(h2, -1.0);
		}
		return diff;
	}

	// create sigma diff histogram
	// bin content = (h1 - h2) / sqrt(h1err^2+h2err^2)

	TH2D* makeRelativeDifferenceHistogram(TH2D *h1, TH2D *h2) const {
		TH2D* reldiff(0);
		if (h1 && h2) {
			reldiff = new TH2D(*h1);
			reldiff->Add(h2, -1.0);
			reldiff->Divide(h2);
		}
		return reldiff;
	}

	TH2D* makeRatioHistogram(TH2D *h1, TH2D *h2) const {
		TH2D* diff(0);
		if (h1 && h2) {
			diff = new TH2D(*h1);
			diff->Divide(h2);
		}
		return diff;
	}

	TH1D* rescaleAxis(TH1D* hist, double scale_factor_x) {
		int nbinsx = hist->GetNbinsX();
		TH1D* hist_rescaled = new TH1D(hist->GetName(), hist->GetTitle(), nbinsx,
				hist->GetXaxis()->GetXmin() * scale_factor_x,
				hist->GetXaxis()->GetXmax() * scale_factor_x);
		for (int ibinx = 1; ibinx <= nbinsx; ibinx++) {
			hist_rescaled->SetBinContent(ibinx, hist->GetBinContent(ibinx));
		}
		return hist_rescaled;
	}
	TH2D* rescaleAxis(TH2D* hist, double scale_factor_x, double scale_factor_y) {
		int nbinsx = hist->GetNbinsX();
		int nbinsy = hist->GetNbinsY();
		TH2D* hist_rescaled = new TH2D(hist->GetName(), hist->GetTitle(), nbinsx,
				hist->GetXaxis()->GetXmin() * scale_factor_x,
				hist->GetXaxis()->GetXmax() * scale_factor_x, nbinsx,
				hist->GetYaxis()->GetXmin() * scale_factor_y,
				hist->GetYaxis()->GetXmax() * scale_factor_y);
		for (int ibinx = 1; ibinx <= nbinsx; ibinx++) {
			for (int ibiny = 1; ibiny <= nbinsy; ibiny++) {
				hist_rescaled->SetBinContent(ibinx, ibiny,
						hist->GetBinContent(ibinx, ibiny));
			}
		}
		return hist_rescaled;
	}

};

class SystematicsAnalyser {
public:
	struct SystematicDependencyGraphBundle {
		TGraphAsymmErrors* mean;
		TGraphAsymmErrors* sigma;
		TGraphAsymmErrors* mean_individual_error; // mean error obtained from individual results
	};

	typedef std::pair<double, double> ValueErrorPair;

private:
	std::map<double, std::vector<double> > dependency_values_map;
	std::map<double, std::vector<double> > dependency_values_error_map;

	struct GausfitPars {
		ValueErrorPair mean;
		ValueErrorPair sigma;
	};

	GausfitPars fitGauss(TH1D *hist) const {
		TF1 gausfit("gauss", "gaus(0)");
		gausfit.SetParameters(10, hist->GetMean(), hist->GetRMS());
		gausfit.SetParNames("A", "#mu", "#sigma");
		hist->Fit(&gausfit, "+");
		GausfitPars gp;
		gp.mean = std::make_pair(gausfit.GetParameter(1), gausfit.GetParError(1));
		gp.sigma = std::make_pair(gausfit.GetParameter(2), gausfit.GetParError(2));
		return gp;
	}

	ValueErrorPair determineMeanAndRMS(const std::vector<double> &values) const {
		double rms(0.0);
		double mean(0.0);
		if (values.size() > 0) {
			for (unsigned int i = 0; i < values.size(); i++) {
				rms = rms + pow(values[i], 2.0);
				mean += values[i];
			}
			rms = sqrt(rms / values.size());
			mean = mean / values.size();
		}
		return std::make_pair(mean, rms);
	}

	std::map<double, TH1D*> createDependencyHistogramMap() const {
		// create histograms
		double sigma_allowance = 4.0;
		std::map<double, TH1D*> dependency_hist_map;

		std::map<double, std::vector<double> >::const_iterator dependency_values;

		for (dependency_values = dependency_values_map.begin();
				dependency_values != dependency_values_map.end(); dependency_values++) {

			// create actual histogram
			std::pair<double, double> hist_bounds = determineMeanAndRMS(
					dependency_values->second);

			dependency_hist_map[dependency_values->first] = new TH1D("dist", "", 100,
					hist_bounds.first - sigma_allowance * hist_bounds.second,
					hist_bounds.first + sigma_allowance * hist_bounds.second);

			// fill data
			const std::vector<double> &value_vec = dependency_values->second;
			for (unsigned int i = 0; i < value_vec.size(); i++) {
				dependency_hist_map[dependency_values->first]->Fill(value_vec[i]);
			}
		}
		return dependency_hist_map;
	}

	std::map<double, GausfitPars> createDependencyGausfitParameterMap(
			const std::map<double, TH1D*> &dependency_hist_map) const {
		std::map<double, GausfitPars> dependency_fit_par_map;

		// gaus fits histograms
		for (std::map<double, TH1D*>::const_iterator it =
				dependency_hist_map.begin(); it != dependency_hist_map.end(); it++) {
			dependency_fit_par_map[it->first] = fitGauss(it->second);
		}

		return dependency_fit_par_map;
	}

	SystematicDependencyGraphBundle createSystematicsGraphBundle(
			const std::map<double, GausfitPars> &dependency_gausfit_par_map) const {
		GraphAndHistogramHelper plot_helper;

		SystematicDependencyGraphBundle graph_bundle;

		std::vector<GraphPoint> mean_vals;
		std::vector<GraphPoint> sigma_vals;
		std::vector<GraphPoint> mean_individual_error_vals;

		GraphPoint temp_graph_point;

		for (std::map<double, GausfitPars>::const_iterator it =
				dependency_gausfit_par_map.begin();
				it != dependency_gausfit_par_map.end(); it++) {

			GausfitPars gausfit_pars = it->second;
			temp_graph_point.x = it->first;

			// mean
			temp_graph_point.y = gausfit_pars.mean.first;
			temp_graph_point.y_err_low = gausfit_pars.mean.second;
			temp_graph_point.y_err_high = gausfit_pars.mean.second;
			mean_vals.push_back(temp_graph_point);
			// sigma
			temp_graph_point.y = gausfit_pars.sigma.first;
			temp_graph_point.y_err_low = gausfit_pars.sigma.second;
			temp_graph_point.y_err_high = gausfit_pars.sigma.second;
			sigma_vals.push_back(temp_graph_point);
		}

		// mean individual error
		for (std::map<double, std::vector<double> >::const_iterator it =
				dependency_values_error_map.begin();
				it != dependency_values_error_map.end(); it++) {
			temp_graph_point.x = it->first;
			ValueErrorPair mean_and_rms = determineMeanAndRMS(it->second);
			temp_graph_point.y = mean_and_rms.first;
			temp_graph_point.y_err_low = mean_and_rms.second;
			temp_graph_point.y_err_high = mean_and_rms.second;
			mean_individual_error_vals.push_back(temp_graph_point);
		}
		graph_bundle.mean_individual_error = plot_helper.makeGraph(
				mean_individual_error_vals);

		return graph_bundle;
	}

public:
	void clear() {
		dependency_values_map.clear();
		dependency_values_error_map.clear();
	}

	void insertDependencyValues(double dependency_value, double measured_value,
			double measured_value_error = 0.0) {
		dependency_values_map[dependency_value].push_back(measured_value);
		if (0.0 != measured_value_error)
			dependency_values_error_map[dependency_value].push_back(
					measured_value_error);
	}

	SystematicDependencyGraphBundle createSystematicsGraphBundle() {
		std::map<double, TH1D*> dependency_hist_map =
				createDependencyHistogramMap();
		std::map<double, GausfitPars> dependency_gausfit_par_map =
				createDependencyGausfitParameterMap(dependency_hist_map);
		return createSystematicsGraphBundle(dependency_gausfit_par_map);
	}

};

} // close namespace

#endif /* ROOTPLOTHELPER_HPP_ */
