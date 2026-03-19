from collections import defaultdict
import json
import os
import traceback

import ipywidgets as widgets
from IPython.display import display
import matplotlib.pylab as plt

from zypotions.widgets import LoginWidget

from session import (
    MutationEffectCalculatorSession,
    USE_LOGICAL_PARENT_STRAIN,
    USE_REFERENCE_STRAIN,
    DEFAULT_CODON_TABLE,
)
import emec.utils.statsmodels_utils as smu
import emec.utils.plot_utils as pu
import emec.utils.translation_utils as tu


class MutationEffectCalculatorInterface:
    """
    This class defines the top-level interface for performing mutation effect analysis for a given
    library of protein mutants.
    """

    # Data inputs keys.
    SEQUENCE_INPUT_MODE = "sequence_input_mode"
    REFERENCE_SEQUENCE_ID = "reference_sequence_id"
    SEQUENCE_FEATURE_TERM = "feature_term"
    SEQUENCE_ANNOTATION_NAME = "annotation_name"
    SEQ_FILENAME_KEY = "sequences_filename"
    SEQUENCE_TRANSLATION_MODE = "sequence_translation_mode"
    SEQUENCE_TRANSLATION_DATA = "sequence_translation_data"
    FIOP_FILENAME_KEY = "fiop_filename"
    FIOP_COL_KEY = "fiop_column"

    # Sequence inputs modes.
    SEQUENCE_INPUT_AUTO = "auto"
    SEQUENCE_INPUT_MANUAL = "manual"

    # Sequence translation modes.
    SEQUENCE_TRANSLATION_STANDARD = "standard"
    SEQUENCE_TRANSLATION_CUSTOM = "custom"

    # Modeling input keys
    PIPELINE_PARAMS_FILENAME_KEY = "pipeline_params_filename"
    DEFAULT_PIPELINE_PARAMS_FILENAME = "pipeline_parameters.json"

    # Reporting inputs
    LIBRARY_TYPE_KEY = "library_type"
    COEFFICIENT_CUTOFF_P_KEY = "coefficient_cutoff_p"
    DEFAULT_COEFFICIENT_CUTOFF_P = 0.001

    # Library types.
    SSL = "Single Site Library"
    MSL = "Combinatorial Library"

    # Default input layout.
    INPUT_LAYOUT = {"width": "650px"}

    def __init__(self):
        """
        Initialize the interface.
        """
        # Setup outputs.
        self.notifications_output = widgets.Output()
        self.results_output = widgets.Output()
        self.save_results_output = widgets.Output()

        # Init user inputs and session.
        self.user_inputs = defaultdict(lambda: None)
        self.session = MutationEffectCalculatorSession()

        # Setup the login widget.
        self.login = LoginWidget()
        self.login.display()

        # Display the interface once the user has logged in.
        self.login.login_button.on_click(
            on_success(self.display, lambda *args: self.login.is_connected)
        )

    def display(self, *args):
        """
        Display the interface.
        """
        title = widgets.HTML(value="<h1>Enzyme Mutation Effects Calculator</h1>")
        display(
            widgets.VBox(
                [
                    title,
                    self.custom_styles,
                    self.form_inputs,
                    self.notifications_output,
                    self.results_output,
                    self.save_results_output,
                ]
            )
        )

    def validate(self):
        """
        Validate input state and display error messages.

        Returns:
            is_valid (bool): True if tab state passes validation.
        """
        validation_errors = []

        # magic number three.
        # TODO(flash): Replace with more nuanced validation.
        if len(self.user_inputs) < 5:
            validation_errors.append(
                "Please provides values for all inputs and try again."
            )

        if validation_errors:
            display_errors(self.notifications_output, validation_errors)

        return len(validation_errors) == 0

    def run(self, widget=None):
        """
        Run the mutation effect analysis.

        Args:
            widget:  optional calling widget if this method is executed as a callback.
        """
        # Clear the notifications_output state.
        self.notifications_output.clear_output()
        self.results_output.clear_output()
        self.save_results_output.clear_output()

        # Validate inputs.
        if not self.validate():
            return

        # Copy user inputs to the inputs directory for posterity.
        with open(os.path.join(self.session.input_dir, "user_inputs.json"), "w") as f:
            json.dump(self.user_inputs, f, indent=True)

        # Load Input data.
        display_notification(self.notifications_output, "Loading user input data...")
        try:
            fiop_fn = self.user_inputs.get(self.FIOP_FILENAME_KEY)
            seq_fn = self.user_inputs.get(self.SEQ_FILENAME_KEY)
            ref_seq_id = self.user_inputs.get(self.REFERENCE_SEQUENCE_ID)
            seq_feature_term = self.user_inputs.get(self.SEQUENCE_FEATURE_TERM)
            seq_annotation_name = self.user_inputs.get(self.SEQUENCE_ANNOTATION_NAME)

            # Set custom codon table to the session for translation.
            tr_mode = self.user_inputs.get(self.SEQUENCE_TRANSLATION_MODE)
            if tr_mode == self.SEQUENCE_TRANSLATION_CUSTOM:
                self.session.codon_table = tu.parse_codon_table(
                    self.user_inputs.get(self.SEQUENCE_TRANSLATION_DATA)
                )
            else:
                self.session.codon_table = DEFAULT_CODON_TABLE

            self.session.load(
                fiop_fn, seq_fn, ref_seq_id, seq_feature_term, seq_annotation_name
            )
        except Exception as er:
            display_notification(
                self.notifications_output,
                f"Error loading data: {str(er)}",
                style="warning",
                show_stacktrace=True,
            )
            return

        # Proceed with analysis.
        display_notification(
            self.notifications_output,
            "Evaluating mutation effects.  This may take a few minutes...",
        )
        try:
            fiop_cn = self.user_inputs.get(self.FIOP_COL_KEY)
            pipeline_params_fn = self.user_inputs.get(self.PIPELINE_PARAMS_FILENAME_KEY)
            self.session.run(fiop_cn, pipeline_params_fn)
        except Exception as er:
            display_notification(
                self.notifications_output,
                f"Error running analysis: {str(er)}",
                style="warning",
                show_stacktrace=True,
            )
            return

        # Display the results.
        display_notification(
            self.notifications_output,
            "Preparing report.  This may take a few minutes...",
        )

        try:
            with self.results_output:
                display(self.results_widget)
        except Exception as er:
            display_notification(
                self.notifications_output,
                f"Error displaying results: {str(er)}",
                style="warning",
                show_stacktrace=True,
            )
            return

    def save(self, button=None):
        """
        Save inputs and results to LIMS as a Dataset.

        Args:
            button: the calling button widget.
        """
        self.save_results_output.clear_output()

        display_notification(
            self.save_results_output, "Saving to Dataset. Please wait..."
        )

        try:
            self.session.save()

            display_notification(
                self.save_results_output,
                widgets.HTML(
                    f"""
                    <p>
                    Data saved to Dataset
                    <a href='{self.session.dataset_link}'
                    target='_blank'>{self.session.dataset_id}</a>
                    </p>"""
                ),
            )
            button.button_style = "success"
        except Exception as er:
            button.button_style = "danger"
            display_notification(
                self.save_results_output,
                f"Error saving results to LIMS: {str(er)}",
                style="danger",
                show_stacktrace=True,
            )

    @property
    def custom_styles(self):
        """
        Returns a css "style" widget, used to customize the interface appearance.
        """
        return widgets.HTML(
            """
            <style>
            .monospace textarea {font-family: monospace;}
            </style>
            """
        )

    @property
    def form_inputs(self):
        """
        Returns the form_inputs widget for the interface.
        """
        upload_prompt = widgets.HTML(
            """
            <h4>Data Inputs</h4>
            <p><i>Upload required input files <a href="./input_files" target="_blank">here</a>.
            Click the "refresh" buttons next to file selectors to see your changes.</i></p>
            """
        )
        return widgets.VBox(
            [
                upload_prompt,
                self.performance_file_input,
                self.sequence_input,
                self.translation_input,
                self.modeling_inputs,
                self.reporting_inputs,
                self.run_button,
            ]
        )

    @property
    def sequence_input(self):
        """
        Returns the sequence_input widget for the interface.

        This widget allows the user to toggle between manual sequence input mode and auto sequence
        input mode.
        """
        # Widget elements.
        button_label = widgets.Label("Sequence Input Mode:")
        button = widgets.RadioButtons(
            options=[self.SEQUENCE_INPUT_AUTO, self.SEQUENCE_INPUT_MANUAL],
            value=self.SEQUENCE_INPUT_AUTO,
            layout=self.INPUT_LAYOUT.copy(),
        )
        auto = self.sequence_input_options_auto
        manual = self.sequence_file_input_manual

        # Auto mode is enabled by default.
        manual.layout.display = "none"
        self.user_inputs[self.SEQUENCE_INPUT_MODE] = self.SEQUENCE_INPUT_AUTO

        # Display controls.
        def toggle_input_display(*args):
            if button.value == self.SEQUENCE_INPUT_MANUAL:
                auto.layout.display = "none"
                manual.layout.display = None
            elif button.value == self.SEQUENCE_INPUT_AUTO:
                auto.layout.display = None
                manual.layout.display = "none"

        # Callbacks.
        button.observe(self.set_state(self.SEQUENCE_INPUT_MODE), names="value")
        button.on_trait_change(toggle_input_display)

        return widgets.VBox([widgets.HBox([button_label, button]), auto, manual])

    @property
    def sequence_input_options_auto(self):
        """
        Returns auto sequence input options.
        """
        ref_id = widgets.Combobox(
            placeholder=(
                "Enter Reference Sequence ID (Strain or DnaComponent), "
                "or select an option from the dropdown."
            ),
            options=[USE_LOGICAL_PARENT_STRAIN, USE_REFERENCE_STRAIN],
            description="Reference Id:",
            layout=self.INPUT_LAYOUT.copy(),
        )
        ref_id.observe(self.set_state(self.REFERENCE_SEQUENCE_ID), names="value")

        ref_at = widgets.Dropdown(
            options=["CDS", "gene"],
            description="Feat. Term:",
            layout=self.INPUT_LAYOUT.copy(),
        )
        ref_at.observe(self.set_state(self.SEQUENCE_FEATURE_TERM, "CDS"), names="value")

        ref_an = widgets.Text(
            placeholder="Enter Annotation Name",
            description="Annot. Name:",
            layout=self.INPUT_LAYOUT.copy(),
        )
        ref_an.observe(self.set_state(self.SEQUENCE_ANNOTATION_NAME), names="value")

        return widgets.VBox([ref_id, ref_at, ref_an])

    @property
    def sequence_file_input_manual(self):
        """
        Returns a sequence file input form widget
        """
        return file_form_input(
            self.session.input_dir,
            "Select Sequences File",
            callback=self.set_state(self.SEQ_FILENAME_KEY),
            file_formats=["fa", "fasta"],
            layout=self.INPUT_LAYOUT.copy(),
        )

    @property
    def translation_input(self):
        """
        Returns the sequence translation mode selector widget.

        This widget allows one to toggle between standard and custom sequence translation
        mode.
        """
        button_label = widgets.Label("Sequence Translation  Mode:")
        button = widgets.RadioButtons(
            options=[
                self.SEQUENCE_TRANSLATION_STANDARD,
                self.SEQUENCE_TRANSLATION_CUSTOM,
            ],
            value=self.SEQUENCE_TRANSLATION_STANDARD,
            layout=self.INPUT_LAYOUT.copy(),
        )

        # Standard mode is enabled by default.
        custom = self.sequence_translation_custom_input
        custom.layout.display = "none"

        # Display controls.
        def toggle_input_display(*args):
            if button.value == self.SEQUENCE_TRANSLATION_STANDARD:
                custom.layout.display = "none"
            elif button.value == self.SEQUENCE_TRANSLATION_CUSTOM:
                custom.layout.display = None

        # Callbacks.
        button.observe(
            self.set_state(
                self.SEQUENCE_TRANSLATION_MODE, self.SEQUENCE_TRANSLATION_STANDARD
            ),
            names="value",
        )
        button.on_trait_change(toggle_input_display)

        return widgets.VBox([widgets.HBox([button_label, button]), custom])

    @property
    def sequence_translation_custom_input(self):
        """
        Returns a custom sequence translation textarea widget.
        """
        translation_data = widgets.Textarea(
            value=tu.STANDARD_CODON_TABLE_DATA,
            description="Codon Table:",
            layout={"width": "800px", "height": "150px"},
        )
        translation_data.add_class("monospace")
        translation_data.observe(
            self.set_state(
                self.SEQUENCE_TRANSLATION_DATA, tu.STANDARD_CODON_TABLE_DATA
            ),
            names="value",
        )

        return translation_data

    @property
    def performance_file_input(self):
        """
        Returns a phenotypes file input form widget.
        """
        ffi = file_form_input(
            self.session.input_dir,
            "Select FIOP File",
            callback=self.set_state(self.FIOP_FILENAME_KEY),
            file_formats=["csv"],
            layout=self.INPUT_LAYOUT.copy(),
        )

        # Attach column identifiers inputs as well.
        # TODO(flash): Populate dynamically from input file.
        cni = widgets.Text(
            placeholder="Enter FIOP column name",
            description="FIOP Col.:",
            layout=self.INPUT_LAYOUT.copy(),
        )
        cni.observe(self.set_state(self.FIOP_COL_KEY), names="value")

        return widgets.VBox([ffi, cni])

    @property
    def modeling_inputs(self):
        """
        Returns modeling parameter widget.
        """
        prompt = widgets.HTML(
            f"""
            <h4>Modeling Parameters</h4>
            <p><i>To change default modeling parameter settings, edit the
            {self.DEFAULT_PIPELINE_PARAMS_FILENAME} file
            <a href="./input_files" target="_blank">here</a>.</i></p>
            """
        )

        ffi = file_form_input(
            self.session.input_dir,
            "Select Modeling Pipeline parameters file",
            callback=self.set_state(
                self.PIPELINE_PARAMS_FILENAME_KEY, self.DEFAULT_PIPELINE_PARAMS_FILENAME
            ),
            file_formats=["json"],
            layout=self.INPUT_LAYOUT.copy(),
        )

        return widgets.VBox([prompt, ffi])

    @property
    def reporting_inputs(self):
        """
        Returns report parameter widget.
        """
        prompt = widgets.HTML("<h4>Reporting Parameters</h4>")

        libd = widgets.Dropdown(
            options=[self.SSL, self.MSL],
            description="Library Type:",
            layout=self.INPUT_LAYOUT.copy(),
        )
        libd.observe(self.set_state(self.LIBRARY_TYPE_KEY, self.SSL), names="value")

        ft = widgets.BoundedFloatText(
            value=self.DEFAULT_COEFFICIENT_CUTOFF_P,
            description="Coef. Cutoff Prob:",
            min=0,
            max=1,
        )
        ft.observe(
            self.set_state(
                self.COEFFICIENT_CUTOFF_P_KEY, self.DEFAULT_COEFFICIENT_CUTOFF_P
            ),
            names="value",
        )

        return widgets.VBox([prompt, libd, ft])

    @property
    def run_button(self):
        """
        Returns a button widget, wired to run the analysis.
        """
        run_button = widgets.Button(
            description="Score mutations",
            layout=self.INPUT_LAYOUT.copy(),
            button_style="primary",
        )
        run_button.on_click(self.run)

        return run_button

    @property
    def results_widget(self) -> widgets.Tab:
        """
        Returns the results widget.
        """
        # Define tabs.
        titles = [
            "Model Summary",
            "Coefficients",
            "Coefficient Interactions",
            "Coefficients Support",
            "Mutation Coverage",
            "Mutation Co-Occurrence",
        ]

        # Define tab contents.
        children = [
            self.model_summary_widget,
            self.model_coefficient_widget,
            widgets.Output(),
            widgets.Output(),
            widgets.Output(),
            widgets.Output(),
        ]
        # Render plots separately to ensure that they make it to the correct output instance.
        with children[2]:
            if self.user_inputs[self.LIBRARY_TYPE_KEY] == self.MSL:
                fn = self.session.output_dir + "/model_interactions.svg"
                p_threshold = self.user_inputs[self.COEFFICIENT_CUTOFF_P_KEY]
                plt.show(
                    pu.plot_model_coefficient_interactions(
                        self.session.results, p_threshold=p_threshold, saveas=fn
                    )
                )
            else:
                msg = f"This tab is disabled. To enable, select the '{self.MSL}' library type."
                display(widgets.HTML(msg))

        # Render plots separately to ensure that they make it to the correct output instance.
        with children[3]:
            fn = self.session.output_dir + "/model_support.svg"
            p_threshold = self.user_inputs[self.COEFFICIENT_CUTOFF_P_KEY]
            plt.show(
                pu.plot_model_coefficient_support(
                    self.session.results, p_threshold=p_threshold, saveas=fn
                )
            )

        # Render plots separately to ensure that they make it to the correct output instance.
        with children[4]:
            fn = self.session.output_dir + "/sequence_coverage.svg"
            plt.show(pu.plot_library_coverage(self.session.alignment, saveas=fn))

        # Render plots separately to ensure that they make it to the correct output instance.
        with children[5]:
            if self.user_inputs[self.LIBRARY_TYPE_KEY] == self.MSL:
                fn = self.session.output_dir + "/mutation_co_occurence.svg"
                plt.show(
                    pu.plot_mutation_co_occurence(self.session.alignment, saveas=fn)
                )
            else:
                msg = f"This tab is disabled. To enable, select the '{self.MSL}' library type."
                display(widgets.HTML(msg))

        # Assemble the tab widget.
        tabs = widgets.Tab(children=children)
        for i, title in enumerate(titles):
            tabs.set_title(i, title)

        # Link to output dir
        download_link = widgets.HTML(
            value=(
                f"""
            <p><i><a href="./output_files/" target="_blank">Download Results To Computer</a>.</i>
                """
            )
        )

        w = widgets.VBox(
            [
                widgets.HTML(value="<h1>Results</h1>"),
                tabs,
                download_link,
                self.save_to_lims_button,
            ]
        )

        return w

    @property
    def model_summary_widget(self) -> widgets.HBox:
        """
        Returns the model summary widget.
        """
        contents = [
            as_widget(
                smu.get_model_summary(
                    self.session.results, self.session.hyperparameters
                ),
                layout=widgets.Layout(width="30%"),
            ),
            widgets.Output(layout=widgets.Layout(width="70%")),
        ]

        with contents[1]:
            plt.show(
                pu.plot_model_fit(
                    self.session.results,
                    saveas=self.session.output_dir + "/model_fit.svg",
                )
            )

        return widgets.HBox(contents)

    @property
    def model_coefficient_widget(self) -> widgets.HBox:
        """
        Returns the model coefficient widget.
        """
        p_threshold = self.user_inputs[self.COEFFICIENT_CUTOFF_P_KEY]

        contents = [
            as_widget(
                smu.get_model_coefficients(
                    self.session.results, p_threshold=p_threshold
                ),
                layout=widgets.Layout(width="43%"),
            ),
            widgets.Output(layout=widgets.Layout(width="57%")),
        ]

        with contents[1]:
            plt.show(
                pu.plot_model_coefficients(
                    self.session.results,
                    p_threshold=p_threshold,
                    saveas=self.session.output_dir + "/model_coefficients.svg",
                )
            )

        return widgets.HBox(contents)

    @property
    def save_to_lims_button(self):
        """
        Returns a button widget, wired to run the analysis.
        """
        save_button = widgets.Button(
            description="Upload Results To LIMS", layout=self.INPUT_LAYOUT.copy()
        )
        save_button.on_click(self.save)

        return save_button

    def set_state(self, state_key: str, default=None):
        """
        This method returns a function that can be used as a widget's observe callback method to
        set update the interface state to the widget's observed value.

        Example:
            select_file = widgets.Select(...)
            select_file.observe(self.set_state("input_file"), names="value")

        Args:
            state_key (str): key in state to store value in.
            default: optional set a default value in case state is not set

        Returns:
            observable callback function.
        """
        if default is not None:
            self.user_inputs[state_key] = default

        def handle_change(change):
            self.user_inputs[state_key] = change.new

        return handle_change


def display_notification(output_widget, message, style="info", show_stacktrace=False):
    """
    Pretty print a message to an output widget.

    Args:
        output_widget (Widget): Output widget to display errors in.
        message (str or widget): Message
        style (str): The message styling, such as "info" or "warning".
        show_stacktrace (bool): Optionally include a stacktrace with the message
    """
    if isinstance(message, str):
        message = message.replace("\n", "<br>")
        children = [widgets.HTML(message)]
    else:
        children = [message]

    if show_stacktrace:
        message = f"<code>{traceback.format_exc()}</code>"
        children.append(widgets.HTML(message))

    val_out = widgets.VBox(children, box_style=style)
    with output_widget:
        display(val_out)


def display_errors(output_widget, errors):
    """
    Helper function for displaying a list of errors in a provided output widget.

    Args:
        output_widget (Widget): Output widget to display errors in.
        errors (List): List of strings of error messages.

    """
    errors_li = [f"<li>{error}</li>" for error in errors]
    error_message = f"""<p>Errors:<ul>{" ".join(errors_li)}</ul></p>"""
    error_message = error_message.replace("\n", "<br/>")
    error_box = widgets.Box([widgets.HTML(error_message)], box_style="warning")
    with output_widget:
        display(error_box)


def file_form_input(
    input_directory, prompt=None, callback=None, file_formats=None, layout=None
):
    """
    Create a drowp-down file form input.

    Args:
        input_directory: the path to the input directory used to populate the dropdown
        prompt: placeholder message and value
        callback: function to call when the widge is observed.
        file_formats: optional list of valid (case-insensitive) file extensions.
        layout: optional layout to apply to the widget

    Returns:
        file dropdown widget.
    """
    prompt = prompt or "Select a File"
    default_option = f"-- {prompt} --"

    # Create dropdown widget.
    input_files_dropdown = widgets.Dropdown(
        options=[default_option], value=default_option, description="Select File:"
    )

    # Define function to update the
    def update_dropdown_options(*args):
        """
        Helper function to update the dropdown widget options based on the input directory contents.
        """
        input_files = [file_name for file_name in os.listdir(input_directory)]
        if file_formats is not None:
            input_files = [
                f
                for f in input_files
                for ext in file_formats
                if f.lower().endswith(ext)
            ]
        input_files_dropdown.options = [default_option] + input_files

    # Update the dropdown options.
    update_dropdown_options()

    # Apply layout if provided.
    if layout is not None:
        input_files_dropdown.layout = layout

    # Attach callback if provided.
    if callback:
        input_files_dropdown.observe(callback, names="value")

    # Define dropdown refresh button
    refresh_button = widgets.Button(
        icon="refresh",
        tooltip="Refresh File List",
        layout=widgets.Layout(width="40px", height="25px"),
    )

    # Call dropdown update on click.
    refresh_button.on_click(update_dropdown_options)

    return widgets.HBox([input_files_dropdown, refresh_button])


def as_widget(d, **kwargs) -> widgets.Output:
    """
    Wrap a displayable object in an Output widget.

    Supported objects include:
        * pandas DataFrame
        * matplotlib Figure
        * matplotlib Axes

    Args:
        d: the displayable object.
        kwargs: optional key-word args to pass to the Output widget.

    Returns:
        an Output widget with display data.
    """
    out = widgets.Output(**kwargs)
    with out:
        display(d)
    return out


def on_success(callback, is_success):
    """
    Execute callback function on success.

    Args:
        callback: the callback function
        is_success: function to check success.
    """

    def f(c):
        if is_success(c):
            callback(c)

    return f
