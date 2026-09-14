use clap::builder::styling::{Ansi256Color, Color, Style, Styles};

const PALE_GOLD: Color = Color::Ansi256(Ansi256Color(222));
const PALE_SKY: Color = Color::Ansi256(Ansi256Color(153));
const PALE_MINT: Color = Color::Ansi256(Ansi256Color(151));
const PALE_APRICOT: Color = Color::Ansi256(Ansi256Color(216));
const SOFT_CORAL: Color = Color::Ansi256(Ansi256Color(210));
const GREY: Color = Color::Ansi256(Ansi256Color(246));

/// A pale headed rosella: the soft yellow of the head on the headings, the wing blue on the flags.
pub const HELP: Styles = Styles::styled()
    .header(Style::new().fg_color(Some(PALE_GOLD)).bold())
    .usage(Style::new().fg_color(Some(PALE_GOLD)).bold())
    .error(Style::new().fg_color(Some(SOFT_CORAL)).bold())
    .literal(Style::new().fg_color(Some(PALE_SKY)))
    .placeholder(Style::new().fg_color(Some(GREY)))
    .valid(Style::new().fg_color(Some(PALE_MINT)))
    .invalid(Style::new().fg_color(Some(PALE_APRICOT)));
