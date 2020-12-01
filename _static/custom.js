/* Adds target=_blank to external links */
$(document).ready(function () {
    $('a[href^="http://"], a[href^="https://"]').not('a[class*=internal]').attr('target', '_blank');
});


/* Formats the default value for each argparse in a code block */
/* <p>Default: <code class="docutils literal notranslate"><span class="pre">10</span></code></p> */
$(document).ready(function () {
    $(".option-list > dd > p").each(function(idx) {
        this.innerHTML = this.innerHTML.replace(/([“”]*)/g, '').replace(/Default: (.+)/g, "Default: <code class=\"docutils literal notranslate\"><span class=\"pre\">$1</span></code>");
    });
});

/* Formats the possible choices for each argparse in a code block */
$(document).ready(function () {
    $(".option-list > dd > p").each(function(idx) {
        let match = this.innerHTML.match(/Possible choices: (.+)/g);
        if (match) {
            let choices = this.innerHTML.replace('Possible choices: ', '').split(', ');
            let formattedChoices = [];
            this.innerHTML = 'Possible choices: '
            for (let i = 0; i < choices.length; i++) {
                formattedChoices.push('<code class="docutils literal notranslate"><span class="pre">' + choices[i] + '</span></code>');
            }
            this.innerHTML = 'Possible choices: ' + formattedChoices.join(', ')
        }
    });
});