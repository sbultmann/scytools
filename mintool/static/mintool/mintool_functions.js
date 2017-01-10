/**
 * Created by bultmann on 1/7/17.
 */
function insert_gRNAs(json) {
    if ($(".gRNA").length != 0) {
        $('.gRNA').remove();
    }
    for (key in json.gRNAs) {
        var gRNA = json.gRNAs[key];
        var table_string = "<div class='gRNA_seq' id='gRNA_seq'>" + gRNA.sequence + "</div>" +
            "<div id='gRNA-desc'>from: " + gRNA.start + " to: " + gRNA.end + "</div>"

        jQuery('<div/>', {
            class: 'gRNA',
            id: key,
            from: gRNA.start,
            to: gRNA.end,
            rel: gRNA.rating,
            strand: gRNA.strand,
            sequence: gRNA.sequence,
            toligo: gRNA.toligo,
            html: '',
        }).appendTo('#gRNA-overview');
        jQuery('<div/>', {
            class: 'gRNA_seq',
            html: gRNA.sequence,
        }).appendTo('#' + key);
        jQuery('<div/>', {
            id: 'gRNA-desc',
            html: "from: " + gRNA.start + " to: " + gRNA.end,
        }).appendTo('#' + key);
    }
    var gRNA_ov = $("#gRNA-overview");
    var gRNA_ov_gRNA = gRNA_ov.children("div");
    gRNA_ov_gRNA.detach().sort(function (a, b) {
        console.log($(a).attr("rel"));
        return $(a).attr('rel') - $(b).attr('rel');
    });

    gRNA_ov.append(gRNA_ov_gRNA);

};

function design_gRNA(sequence, term) {
    console.log("ajax is working!") // sanity check
    $.ajax({
        url: "gRNA-design/", // the endpoint
        type: "GET", // http method
        data: {sequence: sequence, term: term}, // data sent with the post request

        // handle a successful response
        success: function (json) {
            insert_gRNAs(json);
            var data = [];
            for (i = 152, y = 75; i < 852; i = i + 100) {
                data.push([json.sequence.substr(i, 100), y, json.rev_sequence.substr(i, 100), y + 12])
                y = y + 50
            }
            var cords = [bp_to_xy(500)];
            highlight_atg(cords);
            select_gRNA([0,0],1)
            draw_sequence(data);
            draw_ticks();
            $("#primer_f").val(json.primers.left_seq);
            $("#primer_r").val(json.primers.right_seq);
            $("#locus_sequence").val(sequence);
            var primerfwd_pos = [bp_to_xy(json.primers.left_pos[0]), bp_to_xy(json.primers.left_pos[0] + json.primers.left_pos[1])];
            var primerrev_pos = [bp_to_xy(json.primers.right_pos[0] - json.primers.right_pos[1]),bp_to_xy(json.primers.right_pos[0])];
            console.log(primerfwd_pos);
            console.log(primerrev_pos);
            draw_primers(primerfwd_pos, primerrev_pos);
            //console.log(json.sequence); // log the returned json to the console

            console.log("success"); // another sanity check
        },

        // handle a non-successful response
        error: function (xhr, errmsg, err) {
            $('#gRNA-display').html("<div class='alert-box alert radius' data-alert>Oops! We have encountered an error: " + errmsg +
                " <a href='#' class='close'>&times;</a></div>"); // add the error to the dom
            console.log(xhr.status + ": " + xhr.responseText); // provide a bit more info about the error to the console
        }
    });
};

function validateForm() {
  $('.form-field').each(function() {
    if ( $(this).val() === '' ) {
      return false
    }
    else {
      return true;
    }
  });
}


function bp_to_xy(N) {
    N = N - 152;
    var x_bp = N % 100;
    var y_bp = Math.floor(N / 100);
    var y = (y_bp * 50) + 66;
    var x = x_bp * 7;
    return [x, y]
}

function highlight_atg(cords) {
    var container = d3.select("#canvas")
        .attr("height", 500)
        .attr("width", 710);
    var highlight = container.selectAll("#start_stop_codon")
        .data(cords)
        .attr("fill", "red")
        .attr("height", 25)
        .style("opacity", 0.5)
        .attr("id", "start_stop_codon")
        .attr("width", 21)
        .attr("x", function (d) {
            return d[0];
        })
        .attr("y", function (d) {
            return d[1];
        });
    highlight.enter()
        .append("rect")
        .attr("fill", "red")
        .attr("height", 25)
        .style("opacity", 0.5)
        .attr("id", "start_stop_codon")
        .attr("width", 21)
        .attr("x", function (d) {
            return d[0];
        })
        .attr("y", function (d) {
            return d[1];
        });

    highlight.exit().remove();
}

function highlight_gRNA(cords, strand) {
    function split_line(data) {
        var new_data = [];
        new_data.push([data[0][0], [700, data[0][0][1]]]);
        new_data.push([[0, data[0][1][1]], data[0][1]]);
        console.log(new_data);
        return new_data;
    }
    var st = 0;
    var add = 0;
    var num = 1;
    if (strand == 2) {
        add = 12;
        st = 0;
        num = 1;
    }
    else {
        add = 0;
        st = 1;
        num = -1;
    }
    var container = d3.select("#canvas");

    var data_fwd = [cords];
    if (data_fwd[0][1][0] - data_fwd[0][0][0] < 0) {
        data_fwd = split_line(data_fwd);
    }
    var gRNA_highlight = container.selectAll("#gRNA_highlight")
        .data(data_fwd)
        .attr("fill", "grey")
        .attr("height", "12")
        .style("opacity", 0.5)
        .attr("id", "gRNA_highlight")
        .attr("x", function (d) {
            return d[0][0]
        })
        .attr("y", function (d) {
            return d[0][1] + add
        })
        .attr("width", function (d) {
            return d[1][0] - d[0][0]
        });
    gRNA_highlight.enter()
        .append("rect")
        .attr("fill", "grey")
        .attr("height", "12")
        .style("opacity", 0.5)
        .attr("id", "gRNA_highlight")
        .attr("x", function (d) {
            return d[0][0]
        })
        .attr("y", function (d) {
            return d[0][1] + add
        })
        .attr("width", function (d) {
            return d[1][0] - d[0][0]
        });

    gRNA_highlight.exit().remove();

    var cut = container.selectAll("#cut")
        .data(data_fwd)
        .attr("id", "cut")
        .attr("x1", function (d) {
            return d[0][0]+ (49*num) + (st *  (d[1][0] - d[0][0]))
        })
        .attr("x2", function (d) {
            return d[0][0]+ (49*num) + (st *  (d[1][0] - d[0][0]))
        })
        .attr("y1", function (d) {
            return d[0][1] + add - ((12)+ (num * st * 12))
        })
        .attr("y2", function (d) {
            return d[0][1] + add + ((12) + (st * 12))
        })
        .attr("stroke-width", 1)
        .attr("stroke", "black");

    cut.enter()
        .append("line")
        .attr("id", "cut")
        .attr("x1", function (d) {
            return d[0][0]+ 48
        })
        .attr("x2", function (d) {
            return d[0][0]+ 48
        })
        .attr("y1", function (d) {
            return d[0][1] + add - ((12)+ (num * st * 12))
        })
        .attr("y2", function (d) {
            return d[0][1] + add + ((12)+ (st * 12))
        })
        .attr("stroke-width", 1)
        .attr("stroke", "black");
    cut.exit().remove();
}

function select_gRNA(cords, strand) {
    function split_line(data) {
        var new_data = [];
        new_data.push([data[0][0], [700, data[0][0][1]]]);
        new_data.push([[0, data[0][1][1]], data[0][1]]);
        console.log(new_data);
        return new_data;
    }

    if (strand == 2) {
        var add = 12;
    }
    else {
        var add = 0;
    }
    var container = d3.select("#canvas");

    var data_fwd = [cords];
    if (data_fwd[0][1][0] - data_fwd[0][0][0] < 0) {
        data_fwd = split_line(data_fwd);
    }
    var gRNA_select = container.selectAll("#gRNA_select")
        .data(data_fwd)
        .attr("fill", "blue")
        .attr("height", "12")
        .style("opacity", 0.5)
        .attr("id", "gRNA_select")
        .attr("x", function (d) {
            return d[0][0]
        })
        .attr("y", function (d) {
            return d[0][1] + add
        })
        .attr("width", function (d) {
            return d[1][0] - d[0][0]
        });
    gRNA_select.enter()
        .append("rect")
        .attr("fill", "blue")
        .attr("height", "12")
        .style("opacity", 0.5)
        .attr("id", "gRNA_select")
        .attr("x", function (d) {
            return d[0][0]
        })
        .attr("y", function (d) {
            return d[0][1] + add
        })
        .attr("width", function (d) {
            return d[1][0] - d[0][0]
        });

    gRNA_select.exit().remove();
}

function draw_primers(fwd_cords, rev_cords) {
    function split_line_fwd(data) {
        var new_data = [];
        new_data.push([data[0][0], [712, data[0][0][1]]]);
        new_data.push([[0, data[0][1][1]], data[0][1]]);
        console.log(new_data);
        return new_data;
    }
    function split_line_rev(data) {
        var new_data = [];
        new_data.push([data[0][0], [700, data[0][0][1]]]);
        new_data.push([[-12, data[0][1][1]], data[0][1]]);
        console.log(new_data);
        return new_data;
    }

    var container = d3.select("#canvas");

    var data_fwd = [fwd_cords];
    if (data_fwd[0][1][0] - data_fwd[0][0][0] < 0) {
        data_fwd = split_line_fwd(data_fwd);
    }
    var fwd_primer = container.selectAll("#fwd_primer")
        .data(data_fwd)
        .attr("fill", "green")
        .attr("height", "12")
        .style("opacity", 0.5)
        .attr("id", "fwd_primer")
        .attr("x", function (d) {
            return d[0][0]
        })
        .attr("y", function (d) {
            return d[0][1]-12
        })
        .attr("width", function (d) {
            return d[1][0] - d[0][0] -12
        });
    fwd_primer.enter()
        .append("rect")
        .attr("fill", "green")
        .attr("height", "12")
        .style("opacity", 0.5)
        .attr("id", "fwd_primer")
        .attr("x", function (d) {
            return d[0][0]
        })
        .attr("y", function (d) {
            return d[0][1]-12
        })
        .attr("width", function (d) {
            return d[1][0] - d[0][0] -12
        });

    fwd_primer.exit().remove();
    console.log("d", data_fwd[data_fwd.length-1][1]);
    var data_arrow_fwd = [1];
    var lineData = [{"x": data_fwd[data_fwd.length-1][1][0]-12, "y": data_fwd[data_fwd.length-1][1][1]-12},
                    {"x": data_fwd[data_fwd.length-1][1][0]-12, "y": data_fwd[data_fwd.length-1][1][1]-24},
                    {"x": data_fwd[data_fwd.length-1][1][0], "y": data_fwd[data_fwd.length-1][1][1]},
                    {"x": data_fwd[data_fwd.length-1][1][0]-12, "y": data_fwd[data_fwd.length-1][1][1]}];
    console.log("lineData", lineData);
    var lineFunction = d3.line()
        .x(function (d) {
            return d.x;
        })
        .y(function (d) {
            return d.y;
        })
        ;
    var fwd_arrow = container.selectAll("#fwd_arrow")
        .data(data_arrow_fwd)
        .attr("fill", "green")
        .style("opacity", 0.5)
        .attr("id", "fwd_arrow")
        .attr("d", lineFunction(lineData));
    fwd_arrow.enter()
        .append("path")
        .attr("fill", "green")
        .style("opacity", 0.5)
        .attr("id", "fwd_arrow")
        .attr("d", lineFunction(lineData));

    fwd_arrow.exit().remove();


    var data_rev = [rev_cords];
    console.log(data_rev);
    console.log(data_rev);
    if (data_rev[0][1][0] - data_rev[0][0][0] < 0) {
        data_rev = split_line_rev(data_rev);
    }
    var rev_primer = container.selectAll("#rev_primer")
        .data(data_rev)
        .attr("fill", "green")
        .attr("height", "12")
        .style("opacity", 0.5)
        .attr("id", "rev_primer")
        .attr("x", function (d) {
            return d[0][0]+12
        })
        .attr("y", function (d) {
            return d[0][1] + 24
        })
        .attr("width", function (d) {
            return d[1][0] - d[0][0]-12
        });
    rev_primer.enter()
        .append("rect")
        .attr("fill", "green")
        .attr("height", "12")
        .style("opacity", 0.5)
        .attr("id", "rev_primer")
        .attr("x", function (d) {
            return d[0][0]+12
        })
        .attr("y", function (d) {
            return d[0][1] + 24
        })
        .attr("width", function (d) {
            return d[1][0] - d[0][0]-12
        });
    rev_primer.exit().remove();

    var lineData_rev = [{"x": data_rev[0][0][0]+12, "y": data_rev[0][0][1]+36},
                    {"x": data_rev[0][0][0]+12, "y": data_rev[0][0][1]+48},
                    {"x": data_rev[0][0][0], "y": data_rev[0][0][1]+24},
                    {"x": data_rev[0][0][0]+12, "y": data_rev[0][0][1]+24}];
    var rev_arrow = container.selectAll("#rev_arrow")
        .data(data_arrow_fwd)
        .attr("fill", "green")
        .style("opacity", 0.5)
        .attr("id", "rev_arrow")
        .attr("d", lineFunction(lineData_rev));
    rev_arrow.enter()
        .append("path")
        .attr("fill", "green")
        .style("opacity", 0.5)
        .attr("id", "rev_arrow")
        .attr("d", lineFunction(lineData_rev));

    rev_arrow.exit().remove();


}

function draw_ticks() {
    var ticks_data = []
    for (y = 61, n = 0; y < (7 * 50) + 61; y = y + 50) {
        for (i = 133, x = 20; i < 701; i = i + 140) {
            ticks_data.push([i, y, n + x+150]);
            x = x + 20;
        }
        ;
        n = n + 100;
    }
    ;
    var container = d3.select("#canvas");
    var tick = container.selectAll("line")
        .data(ticks_data)
        .enter()
        .append("line")
        .attr("x1", function (d) {
            return d[0] + 3 - 6
        })
        .attr("x2", function (d) {
            return d[0] + 3 - 6
        })
        .attr("y1", function (d) {
            return d[1]
        })
        .attr("y2", function (d) {
            return d[1] + 5
        })
        .attr("stroke-width", 1)
        .attr("stroke", "black");
    var tick_text = container.selectAll("#tick_text")
        .data(ticks_data)
        .enter()
        .append("text")
        .attr("id", "tick_text")
        .attr("font-style", "helvetica")
        .attr("font-size", 8)
        .attr("letter-spacing", "1px")
        .attr("text-anchor", "middle")
        .attr("x", function (d) {
            return d[0] + 3 - 6
        })
        .attr("y", function (d) {
            return d[1] - 3
        })
        .text(function (d) {
            return d[2]
        });
}

function draw_sequence(data) {
    var container = d3.select("#canvas")
        .attr("height", 500)
        .attr("width", 710);
    var fwd = container.selectAll("#fwd")
        .data(data)
        .attr("font-family", "courier new")
        .attr("font-size", 10)
        .attr("letter-spacing", "1px")
        .attr("id", "fwd")
        .attr("x", 0)
        .attr("y", function (d) {
            return d[1]
        })
        .text(function (d) {
            return d[0]
        });
    fwd.enter()
        .append("text")
        .attr("id", "fwd")
        .attr("font-family", "courier new")
        .attr("letter-spacing", "1px")
        .attr("font-size", 10)
        .attr("x", 0)
        .attr("y", function (d) {
            return d[1]
        })
        .text(function (d) {
            return d[0]
        });

    fwd.exit().remove();

    var rev = container.selectAll("#rev")
        .data(data)
        .attr("font-family", "courier new")
        .attr("font-size", 10)
        .attr("letter-spacing", "1px")
        .attr("id", "rev")
        .attr("x", 0)
        .attr("y", function (d) {
            return d[3]
        })
        .text(function (d) {
            return d[2]
        });
    rev.enter()
        .append("text")
        .attr("id", "rev")
        .attr("font-family", "courier new")
        .attr("letter-spacing", "1px")
        .attr("font-size", 10)
        .attr("x", 0)
        .attr("y", function (d) {
            return d[3]
        })
        .text(function (d) {
            return d[2]
        });
    rev.exit().remove();
};