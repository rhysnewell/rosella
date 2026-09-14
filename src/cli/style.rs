use clap::builder::styling::{Ansi256Color, Color, Effects, Style, Styles};

const CRIMSON: Color = Color::Ansi256(Ansi256Color(160));
const STEEL: Color = Color::Ansi256(Ansi256Color(68));
const SAGE: Color = Color::Ansi256(Ansi256Color(107));
const GOLD: Color = Color::Ansi256(Ansi256Color(179));
const GREY: Color = Color::Ansi256(Ansi256Color(245));

/// A crimson rosella: crimson body, the blue of the wing patch on the flags themselves.
pub const HELP: Styles = Styles::styled()
    .header(Style::new().fg_color(Some(CRIMSON)).bold())
    .usage(Style::new().fg_color(Some(CRIMSON)).bold())
    .error(Style::new().fg_color(Some(CRIMSON)).bold())
    .literal(Style::new().fg_color(Some(STEEL)))
    .placeholder(Style::new().fg_color(Some(GREY)))
    .valid(Style::new().fg_color(Some(SAGE)))
    .invalid(Style::new().fg_color(Some(GOLD)));

const _: Effects = Effects::BOLD;
