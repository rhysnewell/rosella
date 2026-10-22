use clap::builder::styling::{Ansi256Color, Color, Style, Styles};

use crate::palette;

const fn tint(colour: u8) -> Option<Color> {
    Some(Color::Ansi256(Ansi256Color(colour)))
}

/// A pale headed rosella: the soft yellow of the head on the headings, the wing blue on the flags.
pub const HELP: Styles = Styles::styled()
    .header(Style::new().fg_color(tint(palette::GOLD)).bold())
    .usage(Style::new().fg_color(tint(palette::GOLD)).bold())
    .error(Style::new().fg_color(tint(palette::CORAL)).bold())
    .literal(Style::new().fg_color(tint(palette::SKY)))
    .placeholder(Style::new().fg_color(tint(palette::GREY)))
    .valid(Style::new().fg_color(tint(palette::MINT)))
    .invalid(Style::new().fg_color(tint(palette::APRICOT)));
